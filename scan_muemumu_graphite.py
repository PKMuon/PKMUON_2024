#!/usr/bin/env python3

import os
import re
import uproot
import numpy as np
import multiprocessing

N_A = 6.022e23
mup_txt_pattern = re.compile(r'^mup_energy_(.*)GeV_Zp_(.*)GeV.txt$')
target_xs = 1e-2 / (2.26 * N_A * 2) * 1e36  # for 2 cm thick graphite, unit pb

mac_template = '''
# ----------------------------------------------
# Settings
# ----------------------------------------------

/control/verbose 2
/gps/verbose 0
/material/verbose 2
/run/verbose 2
/event/verbose 0
/tracking/verbose 0

# ----------------------------------------------
# CLFV Settings
# ----------------------------------------------

/scatter/mupTargetEnToLL {xssf} {mup_txt}
#/run/setCut 0.7 mm

# ----------------------------------------------
# General Particle Source (GPS) Settings
# ----------------------------------------------

/gps/particle mu+
/gps/direction 0.0 0.0 1.0
/gps/totalEnergy {muon_energy} GeV

# ----------------------------------------------
# run
# ----------------------------------------------

/run/initialize
/rlt/SetFileName root_file/mup_{muon_energy}GeV_Zp_{zp_mass}GeV_mue_x{xssf:.0e}.root
/run/printProgress 1000
/run/beamOn 100000
'''[1:]

cnt = 0

def get_xssf(file):
    global cnt
    xses = []
    for rootfile in open(file).read().strip().split():
        print(f'[{cnt}] Reading {rootfile}...')
        cnt += 1
        with uproot.open(rootfile) as file:
            xses.append(file['LHEF']['Event.Weight'].array(library='np').flatten().mean()[0])
    return 10**np.ceil(np.log10(target_xs / np.mean(xses)))

data = []
for file in os.listdir():
    r = mup_txt_pattern.match(file)
    if not r: continue
    muon_energy, zp_mass = r.groups()
    data.append((file, muon_energy, zp_mass))

# Generate mac file.
pool = multiprocessing.Pool()
for (mup_txt, muon_energy, zp_mass), xssf in zip(data, pool.map(get_xssf, [d[0] for d in data])):
    mac_content = mac_template.format(xssf=xssf, mup_txt=mup_txt, muon_energy=muon_energy, zp_mass=zp_mass)
    mac_path = 'build/mup_{muon_energy}GeV_Zp_{zp_mass}GeV_mue_x{xssf:.0e}.mac'.format(
            muon_energy=muon_energy, zp_mass=zp_mass, xssf=xssf)
    print(f'Generating {mac_path}...')
    with open(mac_path, 'w') as mac_file: mac_file.write(mac_content)
