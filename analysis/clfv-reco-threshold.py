#!/usr/bin/env python3

import argparse
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

'''
$ count-events ../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_*.root
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_0.root	997677
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_1.root	997714
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_2.root	997669
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_3.root	997734
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_4.root	997662
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_5.root	997717
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_6.root	997758
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_7.root	997670
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_8.root	997688
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_9.root	997656
$ count-events -s Scatters.Id ../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_*.root
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_0.root	9224
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_1.root	9377
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_2.root	9280
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_3.root	9355
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_4.root	9312
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_5.root	9310
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_6.root	9292
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_7.root	9292
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_8.root	9304
../build/root_file/mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al_9.root	9440
'''
nevent = 997677 + 997714 + 997669 + 997734 + 997662 + 997717 + 997758 + 997670 + 997688 + 997656
nsignal = 9224 + 9377 + 9280 + 9355 + 9312 + 9310 + 9292 + 9292 + 9304 + 9440
nbackground = nevent - nsignal

def savefig(path):
    plt.savefig(path)
    print(f'Plot saved to {path}')

thresholds = np.arange(1, 21) / 10
sig_effs = np.empty_like(thresholds)
bkg_effs = np.empty_like(thresholds)

for i, threshold in enumerate(thresholds):
    rootfile = f'../build/root_file/reco_{threshold}MeV_mup_50.2GeV_Zp_0.25GeV_mumu_x1e+01_al.root'

    # Load input files.
    tree = uproot.concatenate([f'{path}:tree' for path in [rootfile]])
    #print('tree:', *tree.fields, sep='\n  - ')
    for ievent, event in enumerate(tree):
        if ievent >= 2: break
        #print(f'event_{ievent}:', event['Reco.A01'], event['Reco.A02'], event['Reco.A12'], sep='\n  - ')
    
    ## Drop multiple scattering events.
    #tree = tree[ak.num(tree['Scatters.Id'], axis=1) <= 1]
    
    # Split signal and background.
    signal     = tree[tree['MC.IsSignal'] == True ]
    background = tree[tree['MC.IsSignal'] == False]
    for ievent, event in enumerate(signal):
        if ievent >= 10: break
        #print(f'signal_{ievent}:', event['Reco.A01'], event['Reco.A02'], event['Reco.A12'], sep='\n  - ')
    for ievent, event in enumerate(background):
        if ievent >= 10: break
        #print(f'background_{ievent}:', event['Reco.A01'], event['Reco.A02'], event['Reco.A12'], sep='\n  - ')
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

    sig_effs[i] = len(signal) / nsignal
    bkg_effs[i] = len(background) / nbackground

print(thresholds)
print(sig_effs)
print(bkg_effs)
