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
parser.add_argument('-o', metavar='OUTPUT', dest='output', type=str, nargs='?', help='specify output root file')
parser.add_argument('-n', metavar='NEVENT', dest='nevent', type=float, help='specify original event number')
parser.add_argument('-s', metavar='XSSF', dest='xssf', type=float, help='specify global cross section scale factor')
parser.add_argument('-e', metavar='EXPOSURE', dest='exposure', type=float, default=3e13, help='specify total incoming muons')
args = parser.parse_args()
if args.output is None:
    if len(args.input) != 1: raise ValueError('require output name for multiple input names')
    args.output = os.path.join(os.path.dirname(args.input[0]), os.path.basename(args.input[0]))
    args.output = args.output.replace('_*', '').replace('-*', '').replace('*', '').replace('.root', '')
if args.xssf is None:
    if len(args.input) != 1: raise ValueError('require xssf for multiple input names')
    try: args.xssf = float(__import__('re').search(r'[xX]([0-9.eE+-]*)_[^_]*\.root$', args.input[0]).group(1))
    except Exception: raise ValueError('require xssf for unconventional input name')
    print(f'NOTE: XSSF deducted as {args.xssf}')

# Load input files.
tree = uproot.concatenate([f'{path}:tree' for path in args.input])
print('tree:', *tree.fields, sep='\n  - ')
for ievent, event in enumerate(tree):
    if ievent >= 2: break
    print(f'event_{ievent}:', event['Reco.A01'], event['Reco.A02'], event['Reco.A12'], sep='\n  - ')
meta = uproot.concatenate([f'{path}:meta' for path in args.input])
Processes = meta['Meta.Processes'][0]
Processes = Processes.replace('MupTargetEnToLLProcess', r'CLFV-$\mu\mu$')
Processes = Processes.replace('muIoni', r'$\mu$-Ioni')
Processes = Processes.replace('eIoni', r'$e$-Ioni')
Processes = Processes.replace('muPairProd', r'$\mu\mu$-Prod')
Processes = Processes.replace('eBrem', r'$e$-Brem')
Processes = Processes.replace('hadElastic', r'H-Elast')
Processes = Processes.replace('CoulombScat', r'Coulomb')
Processes = Processes.replace('hIoni', r'H-Ioni')
Processes = Processes.replace('pi+Inelastic', r'$\pi^\pm$-Inelast')
Processes = Processes.replace('neutronInelastic', r'$n$-Inelast')
Processes = Processes.replace('annihil', r'Annihil')
Processes = Processes.replace('muBrems', r'$\mu$-Brem')
Processes = Processes.replace('muonNuclear', r'$\mu$-$N$')
Processes = Processes.split('@')
print('processes:', *Processes, sep='\n  - ')

## Drop multiple scattering events.
#tree = tree[ak.num(tree['Scatters.Id'], axis=1) <= 1]

# Compute auxiliary variables.
tree['Reco.A0M'] = np.maximum(tree['Reco.A01'], tree['Reco.A02'])
tree['Reco.AMM'] = np.maximum(tree['Reco.A0M'], tree['Reco.A12'])
p0 = tree['MC.P0']
p1 = tree['MC.P1']
p2 = tree['MC.P2']
up0, up1, up2 = map(lambda p: p / np.linalg.norm(p, axis=-1, keepdims=True), [p0, p1, p2])
tree['MC.A01'] = np.arccos(np.sum(up0 * up1, axis=-1))
tree['MC.A02'] = np.arccos(np.sum(up0 * up2, axis=-1))
tree['MC.A12'] = np.arccos(np.sum(up1 * up2, axis=-1))
tree['MC.A0M'] = np.maximum(tree['MC.A01'], tree['MC.A02'])
tree['MC.AMM'] = np.maximum(tree['MC.A0M'], tree['MC.A12'])

# Split signal and background.
signal     = tree[tree['MC.IsSignal'] == True ]
background = tree[tree['MC.IsSignal'] == False]
for ievent, event in enumerate(signal):
    if ievent >= 10: break
    print(f'signal_{ievent}:', event['Reco.A01'], event['Reco.A02'], event['Reco.A12'], sep='\n  - ')
for ievent, event in enumerate(background):
    if ievent >= 10: break
    print(f'background_{ievent}:', event['Reco.A01'], event['Reco.A02'], event['Reco.A12'], sep='\n  - ')
print('Signal:', ak.num(signal, axis=0))
print('Background:', ak.num(background, axis=0))

def savefig(path):
    plt.savefig(path)
    print(f'Plot saved to {path}')

dof = 0
dof += np.array(tree['Reco.T0']).shape[1] - 2
dof += np.array(tree['Reco.T1']).shape[1] - 2
dof += np.array(tree['Reco.T2']).shape[1] - 2
dof *= 2

# Plot Chi2.
_, bins, _ = plt.hist(signal['Reco.Chi2'], bins=100, range=(0, 10), histtype='step', label='signal')
_, bins, _ = plt.hist(background['Reco.Chi2'], bins=100, range=(0, 10), histtype='step', label='background')
plt.plot(bins, stats.chi2(dof).pdf(bins), label=r'$\chi^2(' + str(dof) + ')$')
plt.xlabel(r'$\chi^2$')
plt.ylabel('Events')
plt.yscale('log')
plt.grid()
plt.legend()
plt.tight_layout()
savefig(args.output + '_Chi2.pdf')
plt.close()

# Apply Chi2 cut.
tree       = tree[tree['Reco.Chi2'] <= dof    ]
signal     = tree[tree['MC.IsSignal'] == True ]
background = tree[tree['MC.IsSignal'] == False]
signal_A0X = np.concatenate([signal['MC.A01'], signal['MC.A02']])

plt.hist(signal['Reco.A01'], bins=100, range=(0, 0.01), histtype='step', label='signal')
plt.hist(background['Reco.A01'], bins=100, range=(0, 0.01), histtype='step', label='background')
plt.hist(signal_A0X, weights=0.5 * np.ones_like(signal_A0X), bins=100, range=(0, 0.01), histtype='step', label='signal-truth')
plt.xlabel(r'<$\vec{p}_0$, $\vec{p}_1$>')
plt.ylabel('Events')
plt.yscale('log')
plt.grid()
plt.legend()
plt.tight_layout()
savefig(args.output + '_A01.pdf')
plt.close()

plt.hist(signal['Reco.A02'], bins=100, range=(0, 0.01), histtype='step', label='signal')
plt.hist(background['Reco.A02'], bins=100, range=(0, 0.01), histtype='step', label='background')
plt.hist(signal_A0X, weights=0.5 * np.ones_like(signal_A0X), bins=100, range=(0, 0.01), histtype='step', label='signal-truth')
plt.xlabel(r'<$\vec{p}_0$, $\vec{p}_2$>')
plt.ylabel('Events')
plt.yscale('log')
plt.grid()
plt.legend()
plt.tight_layout()
savefig(args.output + '_A02.pdf')
plt.close()

plt.hist(signal['Reco.A12'], bins=100, range=(0, 0.01), histtype='step', label='signal')
plt.hist(background['Reco.A12'], bins=100, range=(0, 0.01), histtype='step', label='background')
plt.hist(signal['MC.A12'], bins=100, range=(0, 0.01), histtype='step', label='signal-truth')
plt.xlabel(r'<$\vec{p}_1$, $\vec{p}_2$>')
plt.ylabel('Events')
plt.yscale('log')
plt.grid()
plt.legend()
plt.tight_layout()
savefig(args.output + '_A12.pdf')
plt.close()

plt.hist(signal['Reco.A0M'], bins=100, range=(0, 0.01), histtype='step', label='signal')
plt.hist(background['Reco.A0M'], bins=100, range=(0, 0.01), histtype='step', label='background')
plt.hist(signal['MC.A0M'], bins=100, range=(0, 0.01), histtype='step', label='signal-truth')
plt.xlabel(r'max{<$\vec{p}_0$, $\vec{p}_1$>, <$\vec{p}_0$, $\vec{p}_2$>}')
plt.ylabel('Events')
plt.yscale('log')
plt.grid()
plt.legend()
plt.tight_layout()
savefig(args.output + '_A0M.pdf')
plt.close()

plt.hist(signal['Reco.AMM'], bins=100, range=(0, 0.01), histtype='step', label='signal')
plt.hist(background['Reco.AMM'], bins=100, range=(0, 0.01), histtype='step', label='background')
plt.hist(signal['MC.AMM'], bins=100, range=(0, 0.01), histtype='step', label='signal-truth')
plt.xlabel(r'max{<$\vec{p}_0$, $\vec{p}_1$>, <$\vec{p}_0$, $\vec{p}_2$>, <$\vec{p}_1$, $\vec{p}_2$>}')
plt.ylabel('Events')
plt.yscale('log')
plt.grid()
plt.legend()
plt.tight_layout()
savefig(args.output + '_AMM.pdf')
plt.close()

proc = ak.Array({
    'Processes.Name': Processes,
    'Processes.Contribution': np.mean(signal['Reco.ProcessContributions'], axis=0),
})
proc = proc[proc['Processes.Contribution'] > 0]
proc = sorted([[p['Processes.Contribution'], p['Processes.Name']] for p in proc], reverse=True)
plt.bar(range(len(proc)), [p[0] for p in proc])
plt.xticks(range(len(proc)), [p[1] for p in proc], rotation=45)
plt.xlabel('Process')
plt.ylabel('Mean contribution to energy deposition [MeV]')
plt.yscale('log')
plt.grid(axis='y')
plt.tight_layout()
savefig(args.output + '_Contrib_Sgn.pdf')
plt.close()

proc = ak.Array({
    'Processes.Name': Processes,
    'Processes.Contribution': np.mean(background['Reco.ProcessContributions'], axis=0),
})
proc = proc[proc['Processes.Contribution'] > 0]
proc = sorted([[p['Processes.Contribution'], p['Processes.Name']] for p in proc], reverse=True)
plt.bar(range(len(proc)), [p[0] for p in proc])
plt.xticks(range(len(proc)), [p[1] for p in proc], rotation=45)
plt.xticks(range(len(proc)), [p[1] for p in proc])
plt.xlabel('Process')
plt.ylabel('Mean contribution to energy deposition [MeV]')
plt.yscale('log')
plt.grid(axis='y')
plt.tight_layout()
savefig(args.output + '_Contrib_Bkg.pdf')
plt.close()

if args.nevent is not None:
    signal = signal[signal['Reco.AMM'] < 0.003]
    background = background[background['Reco.AMM'] < 0.003]
    from scipy.stats import chi2
    from scipy.optimize import fsolve
    chi21_95 = chi2(1).ppf(0.95)
    func = lambda s: chi21_95 - 2 * (s + b * np.log(b / (s + b)))
    s = len(signal)
    b = len(background)
    print(f's={s} b={b}')
    if b == 0: b = 1
    s *= args.exposure / args.nevent
    b *= args.exposure / args.nevent
    s_solve = fsolve(func, s)[0]
    diff = func(s_solve)
    assert diff < 1e-3
    ul = args.xssf * s_solve / s
    print(f'ul={ul}')
