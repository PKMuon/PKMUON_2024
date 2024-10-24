#!/usr/bin/env python3

import os
import argparse

# Parse command line arguments.
parser = argparse.ArgumentParser(description='Reconstruct CLFV events.')
parser.add_argument('input', metavar='INPUT', type=str, nargs='+', help='input root files')
parser.add_argument('--output', '-o', type=str, required=True, help='specify output root file')
parser.add_argument('--nevent', '-n', type=float, help='specify the original number of event')
args = parser.parse_args()

import uproot
import awkward as ak
import numba
import numpy as np
from skspatial.objects import Line
import random

# Load input files.
tree = uproot.concatenate([f'{path}:tree' for path in args.input])
meta = uproot.concatenate([f'{path}:meta' for path in args.input])
print('tree:', *tree.fields, sep='\n  - ')
print('meta:', *meta.fields, sep='\n  - ')

# Split signal and background.
signal     = tree[tree['MC.IsSignal'] == True ]
background = tree[tree['MC.IsSignal'] == False]
print('Signal:', ak.num(signal, axis=0))
print('Background:', ak.num(background, axis=0))
tree = background

# Output to ROOT file.
file = uproot.recreate(args.output)
file['tree'] = tree
file['meta'] = {
    'Meta.Processes': np.array([meta[0]['Meta.Processes']]),
    **({'Meta.NEvent': np.array([args.nevent])} if args.nevent is not None else {}),
}
file.close()
print(f'Output written to: {args.output}')
