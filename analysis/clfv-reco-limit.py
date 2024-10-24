#!/usr/bin/env python3

import argparse
import os
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Parse command line arguments.
parser = argparse.ArgumentParser(description='Draw reconstructed CLFV events.')
parser.add_argument('input', metavar='INPUT', type=str, nargs='+', help='input root files')
parser.add_argument('-o', metavar='OUTPUT', dest='output', type=str, required=True, help='specify output root file')
parser.add_argument('-e', metavar='EXPOSURE', dest='exposure', type=float, default=3e13, help='specify total incoming muons')
args = parser.parse_args()

data = { }  # [muon_energy][zp_mass]

for path in args.input:
    # Load input files.
    print(path)
    try:
        muon_energy, zp_mass, xssf = map(float, __import__('re').search(r'mup_([0-9.eE+-]*)GeV_Zp_([0-9.eE+-]*)GeV_mumu_[xX]([0-9.eE+-]*)_[^_]*\.root$', path).groups())
    except Exception:
        print('WARNING: skipping unrecognized path:', path)
        continue
    tree = uproot.concatenate([f'{path}:tree' for path in [path]])
    meta = uproot.concatenate([f'{path}:meta' for path in [path]])
    NEvent = meta['Meta.NEvent'][0]
    tree['Weight'] = args.exposure / NEvent / xssf
    
    # Load supplementary background events.
    import re
    background_path = set(re.sub(r'Zp_.*?GeV_mumu_x[^_]*', 'background_mumu', path) for path in [path])
    if len(background_path) > 2:
        raise RuntimeError(f'Ambiguous background path: {", ".join(sorted(background_path))}')
    if background_path:
        background_path = list(background_path)[0]
        if not os.path.isfile(background_path): background_path = None
    else:
        background_path = None
    print(f'NOTE: supplementary background path: {background_path}')
    if background_path:
        background_tree = uproot.concatenate([f'{background_path}:tree' for _ in range(1)])
        background_meta = uproot.concatenate([f'{background_path}:meta' for _ in range(1)])
        background_NEvent = background_meta['Meta.NEvent'][0]
        background_tree['Weight'] = args.exposure / background_NEvent
        tree = ak.concatenate([tree[tree['MC.IsSignal'] == True], background_tree])
    
    ## Drop multiple scattering events.
    #tree = tree[ak.num(tree['Scatters.Id'], axis=1) <= 1]
    
    # Compute auxiliary variables.
    tree['Reco.A0M'] = np.maximum(tree['Reco.A01'], tree['Reco.A02'])
    tree['Reco.AMM'] = np.maximum(tree['Reco.A0M'], tree['Reco.A12'])
    
    # Split signal and background.
    signal     = tree[tree['MC.IsSignal'] == True ]
    background = tree[tree['MC.IsSignal'] == False]
    print('Signal:', ak.num(signal, axis=0))
    print('Background:', ak.num(background, axis=0))
    
    dof = 0
    dof += np.array(tree['Reco.T0']).shape[1] - 2
    dof += np.array(tree['Reco.T1']).shape[1] - 2
    dof += np.array(tree['Reco.T2']).shape[1] - 2
    dof *= 2
    
    # Apply Chi2 cut.
    tree       = tree[tree['Reco.Chi2'] <= dof    ]
    signal     = tree[tree['MC.IsSignal'] == True ]
    background = tree[tree['MC.IsSignal'] == False]
    
    signal = signal[signal['Reco.AMM'] < 0.003]
    background = background[background['Reco.AMM'] < 0.003]
    from scipy.stats import chi2
    from scipy.optimize import fsolve
    chi21_95 = chi2(1).ppf(0.95)
    func = lambda s: chi21_95 - 2 * (s + b * np.log(b / (s + b)))
    s = np.sum(signal['Weight'])
    b = np.sum(background['Weight'])
    print(f's={s} b={b}')
    if b == 0: b = 1
    s_solve = fsolve(func, s)[0]
    diff = func(s_solve)
    assert diff < 1e-3
    ul = np.sqrt(s_solve / s)
    data[muon_energy] = { **data.get(muon_energy, { }), zp_mass: ul }

def savefig(path):
    plt.savefig(path)
    print(f'Plot saved to {path}')

for muon_energy in data:
    kvp = np.array(sorted(data[muon_energy].items()))
    plt.plot(kvp[:,0], kvp[:,1], label=f'$E_\\mu = {muon_energy:.2f}$ GeV')
plt.xlabel(r'$m_{Z^\prime}$')
plt.ylabel(r'$\lambda_{e\mu}\lambda_{\mu\mu}$ 95% C.L. upper limit')
plt.yscale('log')
plt.grid()
plt.legend()
plt.tight_layout()
savefig(args.output)
plt.close()
