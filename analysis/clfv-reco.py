#!/usr/bin/env python3

import os
import argparse

# Parse command line arguments.
parser = argparse.ArgumentParser(description='Reconstruct CLFV events.')
parser.add_argument('input', metavar='INPUT', type=str, nargs='+', help='input root files')
parser.add_argument('--output', '-o', type=str, nargs='?', help='specify output root file')
parser.add_argument('--electron-veto-rate', '-e', type=float, default=0.0, help='specify electron veto rate')
parser.add_argument('--nevent', '-n', type=float, help='specify the original number of event')
args = parser.parse_args()
if args.output is None:
    if len(args.input) != 1: raise ValueError('require output name for multiple input names')
    args.output = os.path.join(os.path.dirname(args.input[0]), 'reco_' + os.path.basename(args.input[0]))
    args.output = args.output.replace('_*', '').replace('-*', '').replace('*', '')

import uproot
import awkward as ak
import numba
import numpy as np
from skspatial.objects import Line
import random

# Load input files.
tree = uproot.concatenate([f'{path}:tree' for path in args.input])
params = uproot.concatenate([f'{path}:params' for path in args.input])
print('tree:', *tree.fields, sep='\n  - ')
print('params:', *params.fields, sep='\n  - ')
for ievent, event in enumerate(tree):
    if ievent >= 2: break
    print(f'event_{ievent}:', event['Edeps.Id'], event['Edeps.Pid'], event['Edeps.Value'], sep='\n  - ')

# Read parameter setting.
def reduce_param(param):  # [NRun, NField]
    param = param
    assert ak.all(param == param[:1])
    return param[0]
LayerZ = reduce_param(params['Params.LayerZ'])[0]
CellX = reduce_param(params['Params.CellX'])[0]
CellY = reduce_param(params['Params.CellY'])[0]
HalfNCellX = reduce_param(params['Params.HalfNCellX'])[0]
HalfNCellY = reduce_param(params['Params.HalfNCellY'])[0]
Processes = reduce_param(params['Processes.Name'])
print('LayerZ:', LayerZ)
print('CellX:', CellX)
print('CellY:', CellY)
print('HalfNCellX:', HalfNCellX)
print('HalfNCellY:', HalfNCellY)
print('Processes:', Processes)
NLayer = len(LayerZ)
NCellX = 2 * HalfNCellX
NCellY = 2 * HalfNCellY
XOffset = -HalfNCellX * CellX
YOffset = -HalfNCellY * CellY
VarInv = np.array([
    12 / CellX**2,
    12 / CellY**2,
    1.0,  # not inf to avoid (0.0 * inf)
])
NProcess = len(Processes)

# Simulate PID specific detector response.
electron_veto_rate = args.electron_veto_rate
@numba.jit
def simulate_pid_response(tree):  # [NOTE] Edeps.Id must be sorted.
    Edeps_Id, Edeps_Value = [ ], [ ]
    for ievent, event in enumerate(tree):
        edeps_id, edeps_value = [ ], [ ]
        last_id, last_edep = None, 0.0
        for id, pid, edep in zip(event['Edeps.Id'], event['Edeps.Pid'], event['Edeps.Value']):
            if abs(id) == 11 and random.random() < electron_veto_rate: continue
            if last_id is not None and last_id != id: edeps_id.append(last_id); edeps_value.append(last_edep); last_edep = 0.0
            last_id = id; last_edep += edep
        if last_id is not None: edeps_id.append(last_id); edeps_value.append(last_edep)
        Edeps_Id.append(edeps_id)
        Edeps_Value.append(edeps_value)
    return Edeps_Id, Edeps_Value
tree['Edeps.ProcessValue'] = tree['Edeps.Value']
Edeps_Ids = []
Edeps_Values = []
for i in range(0, len(tree), 10000):
    Edeps_Id, Edeps_Value = simulate_pid_response(tree[i : i + 10000])
    Edeps_Ids.append(ak.Array(Edeps_Id))
    Edeps_Values.append(ak.Array(Edeps_Value))
    del Edeps_Id, Edeps_Value
tree['Edeps.Id'] = ak.concatenate(Edeps_Ids)
tree['Edeps.Value'] = ak.concatenate(Edeps_Values)
del Edeps_Ids, Edeps_Values

# Decode energy deposit positions.
tree['Edeps.Layer'] = tree['Edeps.Id'] // (NCellX * NCellY)
@numba.jit
def decode_xyz(tree):
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
Xs = []
Ys = []
Zs = []
for i in range(0, len(tree), 10000):
    subtree = tree[i : i + 10000]
    X, Y, Z = decode_xyz(subtree)
    Xs.append(ak.Array(X))
    Ys.append(ak.Array(Y))
    Zs.append(ak.Array(Z))
    del X, Y, Z
tree['Edeps.X'] = ak.concatenate(Xs)
tree['Edeps.Y'] = ak.concatenate(Ys)
tree['Edeps.Z'] = ak.concatenate(Zs)
del Xs, Ys, Zs
#for ievent, event in enumerate(tree):
#    if ievent >= 10: break
#    print(f'event_{ievent}:', event['Edeps.Layer'], event['Edeps.X'], event['Edeps.Y'], event['Edeps.Value'], sep='\n  - ')

# Simulate readout structure.
assert NLayer % 2 == 0
NLayer //= 2
Edeps = tree[[field for field in tree.fields if field.startswith('Edeps.')]]
Edeps = Edeps[Edeps['Edeps.Value'] >= 1]  # Require minimum energy deposition 1 MeV.
XEdeps = Edeps[Edeps['Edeps.Layer'] % 2 == 0]  # Require consistent x-y waveforms.
YEdeps = Edeps[Edeps['Edeps.Layer'] % 2 == 1]
mask_0 = ak.num(XEdeps['Edeps.Id'], axis=1) == ak.num(YEdeps['Edeps.Id'], axis=1)
XEdeps, YEdeps = XEdeps[mask_0], YEdeps[mask_0]
mask_1 = ak.all(XEdeps['Edeps.Id'] % (NCellX * NCellY) == YEdeps['Edeps.Id'] % (NCellX * NCellY), axis=1)
XEdeps, YEdeps = XEdeps[mask_1], YEdeps[mask_1]
Edeps = XEdeps
Edeps['Edeps.Layer'] = Edeps['Edeps.Layer'] // 2
Edeps['Edeps.Z'] = (XEdeps['Edeps.Z'] + YEdeps['Edeps.Z']) / 2
Edeps['Edeps.Value'] = XEdeps['Edeps.Value'] + YEdeps['Edeps.Value']
tree = tree[mask_0]
tree = tree[mask_1]
for field in Edeps.fields: tree[field] = Edeps[field]
del Edeps, XEdeps, YEdeps, mask_0, mask_1
#for ievent, event in enumerate(tree):
#    if ievent >= 10: break
#    print(f'event_{ievent}:', event['Edeps.Layer'], event['Edeps.X'], event['Edeps.Y'], event['Edeps.Value'], sep='\n  - ')

# Select events with 1 hit in the first half and 2 hits in the second half.
assert NLayer % 2 == 0
hits = sorted([l for l in range(NLayer)] + [l for l in range(NLayer // 2, NLayer)])
tree = tree[ak.num(tree['Edeps.Layer'], axis=1) == len(hits)]
tree = tree[ak.all(tree['Edeps.Layer'] == [hits], axis=1)]

# Reconstruct positions and angles.
assert NLayer >= 4
r = np.array([tree['Edeps.X'].to_numpy(), tree['Edeps.Y'].to_numpy(), tree['Edeps.Z'].to_numpy()])  # [3, NEvent, len(hits)]
r = np.transpose(r, [1, 2, 0])  # [NEvent, len(hits), 3]
r0 = r[:, 0]
r1 = r[:, NLayer // 2]
r2 = r[:, NLayer // 2 + 1]
ra = r[:, len(hits) - 2]
rb = r[:, len(hits) - 1]
d0 = r[:, NLayer // 2 - 1] - r0
d1_a = ra - r1
d2_a = rb - r2
d1_b = rb - r1
d2_b = ra - r2
dist_a = np.linalg.norm(d1_a, axis=-1, keepdims=True) + np.linalg.norm(d2_a, axis=-1, keepdims=True)
dist_b = np.linalg.norm(d1_b, axis=-1, keepdims=True) + np.linalg.norm(d2_b, axis=-1, keepdims=True)
d1 = np.where(dist_a < dist_b, d1_a, d1_b)
d2 = np.where(dist_a < dist_b, d2_a, d2_b)

# Minimal Chi2 matching.
# Needs: r, r0, r1, r2, d0, d1, d2
@numba.jit
def min_chi2_match():
    Chi2 = [ ]
    T0, T1, T2 = [ ], [ ], [ ]
    for ievent, event in enumerate(tree):
        chi2 = 0.0
        t0, t1, t2 = [ ], [ ], [ ]

        # The first half.
        for l in range(NLayer // 2):
            ra = r[ievent, l]
            r0_a = r0[ievent] + d0[ievent] * (ra[2] - r0[ievent, 2]) / d0[ievent, 2]
            chi2 += np.sum((ra - r0_a)**2 * VarInv); t0.append(ra)

        # The second half.
        for l in range(NLayer // 2, NLayer):
            ra = r[ievent, l * 2 - NLayer // 2]
            rb = r[ievent, l * 2 - NLayer // 2 + 1]
            r1_a = r1[ievent] + d1[ievent] * (ra[2] - r1[ievent, 2]) / d1[ievent, 2]
            r2_a = r2[ievent] + d2[ievent] * (rb[2] - r2[ievent, 2]) / d2[ievent, 2]
            r1_b = r1[ievent] + d1[ievent] * (rb[2] - r1[ievent, 2]) / d1[ievent, 2]
            r2_b = r2[ievent] + d2[ievent] * (ra[2] - r2[ievent, 2]) / d2[ievent, 2]
            chi2_a = np.sum((ra - r1_a)**2 * VarInv) + np.sum((rb - r2_a)**2 * VarInv)
            chi2_b = np.sum((rb - r1_b)**2 * VarInv) + np.sum((ra - r2_b)**2 * VarInv)
            if chi2_a <= chi2_b:
                chi2 += chi2_a; t1.append(ra); t2.append(rb)
            else:
                chi2 += chi2_b; t1.append(rb); t2.append(ra)

        Chi2.append(chi2)
        T0.append(t0); T1.append(t1); T2.append(t2)
    return Chi2, T0, T1, T2
Chi2, T0, T1, T2 = map(np.array, min_chi2_match())
tree['Reco.T0'] = T0
tree['Reco.T1'] = T1
tree['Reco.T2'] = T2
#for ievent, event in enumerate(tree):
#    if ievent >= 10: break
#    print(f'event_{ievent}:', event['Reco.T0'], event['Reco.T1'], event['Reco.T2'], sep='\n  - ')
del r, r0, r1, r2, ra, rb, d0, d1_a, d2_a, d1_b, d2_b, dist_a, dist_b, d1, d2

# 3D line fit.
def line_fit(t):  # [NEvent, NLayer // 2, 3]
    lines = list(map(Line.best_fit, t))
    points = np.array([line.point for line in lines])
    directions = np.array([line.direction for line in lines])
    return points, directions
P0, D0 = line_fit(tree['Reco.T0'])
P1, D1 = line_fit(tree['Reco.T1'])
P2, D2 = line_fit(tree['Reco.T2'])
tree['Reco.P0'], tree['Reco.D0'] = P0, D0
tree['Reco.P1'], tree['Reco.D1'] = P1, D1
tree['Reco.P2'], tree['Reco.D2'] = P2, D2
tree['Reco.A01'] = np.arccos(np.minimum(1.0, np.sum(np.abs(tree['Reco.D0'] * tree['Reco.D1']), axis=-1)))
tree['Reco.A02'] = np.arccos(np.minimum(1.0, np.sum(np.abs(tree['Reco.D0'] * tree['Reco.D2']), axis=-1)))
tree['Reco.A12'] = np.arccos(np.minimum(1.0, np.sum(np.abs(tree['Reco.D1'] * tree['Reco.D2']), axis=-1)))

# Recompute Chi2.
P0, D0, P1, D1, P2, D2 = map(lambda x: x.reshape(-1, 1, 3), [P0, D0, P1, D1, P2, D2])
Chi2 = 0
Chi2 += np.sum((P0 + D0 * (T0[:,:,2:] - P0[:,:,2:]) / D0[:,:,2:] - T0)**2 * VarInv.reshape(1, 1, 3), axis=(1, 2))
Chi2 += np.sum((P1 + D1 * (T1[:,:,2:] - P1[:,:,2:]) / D1[:,:,2:] - T1)**2 * VarInv.reshape(1, 1, 3), axis=(1, 2))
Chi2 += np.sum((P2 + D2 * (T2[:,:,2:] - P2[:,:,2:]) / D2[:,:,2:] - T2)**2 * VarInv.reshape(1, 1, 3), axis=(1, 2))
tree['Reco.Chi2'] = Chi2

# Compute process_contributions.
Edeps = tree[['Edeps.Process', 'Edeps.ProcessValue']]
ProcessContributions = np.empty((len(tree), NProcess + 1))
for process in range(-1, NProcess):
    edeps = Edeps[Edeps['Edeps.Process'] == process]
    # [TODO] Filter pid.
    ProcessContributions[:, process] = ak.sum(edeps['Edeps.ProcessValue'], axis=1)
tree['Reco.ProcessContributions'] = ProcessContributions
Processes = ak.Array([*Processes, 'Primary'])
NProcess += 1

# Drop multiple scattering events.
tree = tree[ak.num(tree['Scatters.Id'], axis=1) <= 1]

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
    if ievent >= 2: break
    print(f'signal_{ievent}:', event['Edeps.Layer'], event['Edeps.X'], event['Edeps.Y'],
          event['Edeps.Value'], event['Reco.A12'], event['Reco.Chi2'], sep='\n  - ')
for ievent, event in enumerate(background):
    if ievent >= 2: break
    print(f'background_{ievent}:', event['Edeps.Layer'], event['Edeps.X'], event['Edeps.Y'],
          event['Edeps.Value'], event['Reco.A12'], event['Reco.Chi2'], sep='\n  - ')
print('Signal:', ak.num(signal, axis=0))
print('Background:', ak.num(background, axis=0))

# Output to ROOT file.
tree = tree[[
    field for field in tree.fields
    if field.startswith('Reco.') or field.startswith('MC.')
]]
file = uproot.recreate(args.output)
file['tree'] = tree
file['meta'] = {
    'Meta.Processes': np.array(['@'.join(Processes)]),
    **({'Meta.NEvent': np.array([args.nevent])} if args.nevent is not None else {}),
}
file.close()
print(f'Output written to: {args.output}')
