#!/usr/bin/env python3

import argparse
import os
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

# Parse command line arguments.
parser = argparse.ArgumentParser(description='Draw reconstructed CLFV events.')
parser.add_argument('input', metavar='INPUT', type=str, nargs='+', help='input root files')
parser.add_argument('-o', metavar='OUTPUT', dest='output', type=str, nargs='?', help='specify output root file')
args = parser.parse_args()
if args.output is None:
    if len(args.input) != 1: raise ValueError('require output name for multiple input names')
    args.output = os.path.join(os.path.dirname(args.input[0]), os.path.basename(args.input[0]))
    args.output = args.output.replace('_*', '').replace('-*', '').replace('*', '').replace('.root', '')

# Load input files.
tree = uproot.concatenate([f'{path}:tree' for path in args.input])
print('tree:', *tree.fields, sep='\n  - ')
for ievent, event in enumerate(tree):
    if ievent >= 2: break
    print(f'event_{ievent}:', event['Reco.A01'], event['Reco.A02'], event['Reco.A12'], sep='\n  - ')

## Drop multiple scattering events.
#tree = tree[ak.num(tree['Scatters.Id'], axis=1) <= 1]

# Compute auxiliary variables.
tree['Reco.A0M'] = np.maximum(tree['Reco.A01'], tree['Reco.A02'])
tree['Reco.AMM'] = np.maximum(tree['Reco.A0M'], tree['Reco.A12'])

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

plt.hist(signal['Reco.A01'], bins=10, range=(0, 0.01), histtype='step', label='signal')
plt.hist(background['Reco.A01'], bins=10, range=(0, 0.01), histtype='step', label='background')
plt.xlabel(r'<$\vec{p}_0$, $\vec{p}_1$>')
plt.ylabel('Events')
plt.legend()
plt.tight_layout()
savefig(args.output + '_A01.pdf')
plt.close()

plt.hist(signal['Reco.A02'], bins=10, range=(0, 0.01), histtype='step', label='signal')
plt.hist(background['Reco.A02'], bins=10, range=(0, 0.01), histtype='step', label='background')
plt.xlabel(r'<$\vec{p}_0$, $\vec{p}_2$>')
plt.ylabel('Events')
plt.legend()
plt.tight_layout()
savefig(args.output + '_A02.pdf')
plt.close()

plt.hist(signal['Reco.A12'], bins=10, range=(0, 0.01), histtype='step', label='signal')
plt.hist(background['Reco.A12'], bins=10, range=(0, 0.01), histtype='step', label='background')
plt.xlabel(r'<$\vec{p}_1$, $\vec{p}_2$>')
plt.ylabel('Events')
plt.legend()
plt.tight_layout()
savefig(args.output + '_A12.pdf')
plt.close()

plt.hist(signal['Reco.A0M'], bins=10, range=(0, 0.01), histtype='step', label='signal')
plt.hist(background['Reco.A0M'], bins=10, range=(0, 0.01), histtype='step', label='background')
plt.xlabel(r'min{<$\vec{p}_0$, $\vec{p}_1$>, <$\vec{p}_0$, $\vec{p}_2$>}')
plt.ylabel('Events')
plt.legend()
plt.tight_layout()
savefig(args.output + '_A0M.pdf')
plt.close()

plt.hist(signal['Reco.AMM'], bins=10, range=(0, 0.01), histtype='step', label='signal')
plt.hist(background['Reco.AMM'], bins=10, range=(0, 0.01), histtype='step', label='background')
plt.xlabel(r'min{<$\vec{p}_0$, $\vec{p}_1$>, <$\vec{p}_0$, $\vec{p}_2$>, <$\vec{p}_1$, $\vec{p}_2$>}')
plt.ylabel('Events')
plt.legend()
plt.tight_layout()
savefig(args.output + '_AMM.pdf')
plt.close()
