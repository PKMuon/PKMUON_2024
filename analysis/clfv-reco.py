#!/usr/bin/env python3

import argparse
import os
import uproot
import awkward as ak
import numba
import numpy as np

# Parse command line arguments.
parser = argparse.ArgumentParser(description='Reconstruct CLFV events.')
parser.add_argument('input', metavar='INPUT', type=str, nargs='+', help='input root files')
parser.add_argument('-o', metavar='OUTPUT', dest='output', type=str, nargs='?', help='specify output root file')
args = parser.parse_args()
if args.output is None:
    if len(args.input) != 1: raise ValueError('require output name for multiple input names')
    args.output = os.path.join(os.path.dirname(args.input[0]), 'reco_' + os.path.basename(args.input[0]))
    args.output = args.output.replace('_*', '').replace('-*', '').replace('*', '')

# Load input files.
tree = uproot.concatenate([f'{path}:tree' for path in args.input])
params = uproot.concatenate([f'{path}:params' for path in args.input])
print('tree:', *tree.fields, sep='\n  - ')
print('params:', *params.fields, sep='\n  - ')
for ievent, event in enumerate(tree):
    if ievent >= 2: break
    print(f'event_{ievent}:', event['Edeps.Id'], event['Edeps.Pid'], event['Edeps.Value'], sep='\n  - ')

# Read parameter setting.
def reduce_param(param):  # [NRun, 1]
    param = param
    assert ak.all(param == param[:1])
    return param[0,0]
LayerZ = reduce_param(params['Params.LayerZ'])
CellX = reduce_param(params['Params.CellX'])
CellY = reduce_param(params['Params.CellY'])
HalfNCellX = reduce_param(params['Params.HalfNCellX'])
HalfNCellY = reduce_param(params['Params.HalfNCellY'])
print('LayerZ:', LayerZ)
print('CellX:', CellX)
print('CellY:', CellY)
print('HalfNCellX:', HalfNCellX)
print('HalfNCellY:', HalfNCellY)
NLayer = len(LayerZ)
NCellX = 2 * HalfNCellX
NCellY = 2 * HalfNCellY
XOffset = -HalfNCellX * CellX
YOffset = -HalfNCellY * CellY

# Simulate PID specific detector response.
@numba.jit
def simulate_pid_response():  # [NOTE] Edeps.Id must be sorted.
    Edeps_Id, Edeps_Value = [ ], [ ]
    for ievent, event in enumerate(tree):
        edeps_id, edeps_value = [ ], [ ]
        last_id, last_edep = None, 0.0
        for id, pid, edep in zip(event['Edeps.Id'], event['Edeps.Pid'], event['Edeps.Value']):
            # [TODO] Filter pid.
            if last_id is not None and last_id != id: edeps_id.append(last_id); edeps_value.append(last_edep); last_edep = 0.0
            last_id = id; last_edep += edep
        if last_id is not None: edeps_id.append(last_id); edeps_value.append(last_edep)
        Edeps_Id.append(edeps_id)
        Edeps_Value.append(edeps_value)
    return Edeps_Id, Edeps_Value
Edeps_Id, Edeps_Value = simulate_pid_response()
tree['Edeps.Id'] = Edeps_Id
tree['Edeps.Pid'] = ak.zeros_like(Edeps_Id)
tree['Edeps.Value'] = Edeps_Value

# Decode energy deposit positions.
tree['Edeps.Layer'] = tree['Edeps.Id'] // (NCellX * NCellY)
@numba.jit
def decode_xyz():
    X, Y, Z = [ ], [ ], [ ]
    for ievent, event in enumerate(tree):
        x, y, z = [ ], [ ], [ ]
        for id in event['Edeps.Id']:
            z.append(LayerZ[id // (NCellX * NCellY)])
            id %= NCellX * NCellY
            x.append(XOffset + (id // NCellY + 0.5) * CellX)
            y.append(YOffset + (id % NCellY + 0.5) * CellY)
        X.append(x); Y.append(y); Z.append(z)
    return X, Y, Z
X, Y, Z = decode_xyz()
tree['Edeps.X'] = X; tree['Edeps.Y'] = Y; tree['Edeps.Z'] = Z

# Simulate readout structure.
assert NLayer % 2 == 0
NLayer //= 2
Edeps = tree[[field for field in tree.fields if field.startswith('Edeps.')]]
Edeps = Edeps[Edeps['Edeps.Layer'] % 2 == 0]  # [XXX] We drop half of the layers temporarily.
Edeps = Edeps[Edeps['Edeps.Value'] >= 1]  # Require minimum energy deposition 1 MeV.
Edeps['Edeps.Layer'] = Edeps['Edeps.Layer'] // 2
for field in Edeps.fields: tree[field] = Edeps[field]

# Select events with hits pattern [0, 1, 2, 2, 3, 3].
tree = tree[ak.num(tree['Edeps.Layer'], axis=1) == 6]
tree = tree[ak.all(tree['Edeps.Layer'] == [[0, 1, 2, 2, 3, 3]], axis=1)]

# Reconstruct positions and angles.
r = np.array([tree['Edeps.X'].to_numpy(), tree['Edeps.Y'].to_numpy(), tree['Edeps.Z'].to_numpy()])  # [3, NEvent, 6]
r = np.transpose(r, [1, 2, 0])  # [NEvent, 6, 3]
r0 = r[:,0]
r1 = r[:,2]
r2 = r[:,3]
dr0 = r[:,1] - r[:,0]
dr1_a = r[:,4] - r[:,2]
dr2_a = r[:,5] - r[:,3]
dr1_b = r[:,5] - r[:,2]
dr2_b = r[:,4] - r[:,3]
dist_a = np.linalg.norm(dr1_a, axis=-1, keepdims=True) + np.linalg.norm(dr2_a, axis=-1, keepdims=True)
dist_b = np.linalg.norm(dr1_b, axis=-1, keepdims=True) + np.linalg.norm(dr2_b, axis=-1, keepdims=True)
dr1 = np.where(dist_a < dist_b, dr1_a, dr1_b)
dr2 = np.where(dist_a < dist_b, dr2_a, dr2_b)
dr0, dr1, dr2 = map(lambda r: r / np.linalg.norm(r, axis=-1, keepdims=True), [dr0, dr1, dr2])
a01 = np.arccos(np.sum(dr0 * dr1, axis=-1))
a02 = np.arccos(np.sum(dr0 * dr2, axis=-1))
a12 = np.arccos(np.sum(dr1 * dr2, axis=-1))
tree['Reco.R0'] = r0
tree['Reco.R1'] = r1
tree['Reco.R2'] = r2
tree['Reco.D0'] = dr0
tree['Reco.D1'] = dr1
tree['Reco.D2'] = dr2
tree['Reco.A01'] = a01
tree['Reco.A02'] = a02
tree['Reco.A12'] = a12

## Drop multiple scattering events.
#tree = tree[ak.num(tree['Scatters.Id'], axis=1) <= 1]

# Store MC truth.
tree['MC.IsSignal'] = np.where(ak.num(tree['Scatters.Id'], axis=1) == 1, True, False)
r = [tree['Scatters.X'], tree['Scatters.Y'], tree['Scatters.Z']]  # [3, NEvent, PAD]
r = ak.fill_none(ak.pad_none(r, 1, axis=2), 0, axis=2)
p = [tree['Scatters.Px[3]'], tree['Scatters.Py[3]'], tree['Scatters.Pz[3]']]  # [3, NEvent, PAD, 3]
p = ak.fill_none(ak.pad_none(p, 1, axis=2), [0, 0, 0], axis=2)
r = np.transpose(r, [1, 2, 0])  # [NEvent, 1, 3]
p = np.transpose(p, [1, 2, 3, 0])  # [NEvent, 1, 3, 3]
tree['MC.R'] = r[:,0]
tree['MC.P0'] = p[:,0,0]
tree['MC.P1'] = p[:,0,1]
tree['MC.P2'] = p[:,0,2]

# Split signal and background.
signal     = tree[tree['MC.IsSignal'] == True ]
background = tree[tree['MC.IsSignal'] == False]
for ievent, event in enumerate(signal):
    if ievent >= 10: break
    print(f'signal_{ievent}:', event['Edeps.Layer'], event['Edeps.X'], event['Edeps.Y'], event['Edeps.Pid'],
          event['Edeps.Value'], event['Reco.A12'], sep='\n  - ')
for ievent, event in enumerate(background):
    if ievent >= 10: break
    print(f'background_{ievent}:', event['Edeps.Layer'], event['Edeps.X'], event['Edeps.Y'], event['Edeps.Pid'],
          event['Edeps.Value'], event['Reco.A12'], sep='\n  - ')
print('Signal:', ak.num(signal, axis=0))
print('Background:', ak.num(background, axis=0))

# Output to ROOT file.
tree = tree[[field for field in tree.fields if field.startswith('Reco.') or field.startswith('MC.')]]
file = uproot.recreate(args.output)
file['tree'] = tree
file.close()
print(f'Output written to: {args.output}')
