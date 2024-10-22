#!/usr/bin/env python3

import re
import glob
import uproot
import numpy as np
from scipy.constants import N_A
from scipy.interpolate import make_interp_spline

ki_points = glob.glob('clfv/data/muemumu_ki_*.root')
ki_point_pattern = re.compile(r'muemumu_ki_([^_/]*)\.root$')
muon_energies = np.array([33.6, 50.2, 77.2])  # unit GeV
target_xs = 1e-3 / (2.7 * N_A * 8) * 1e36  # for 8 cm thick Al, unit pb

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

/scatter/mupTargetEnToLL 13 {ki_point} {xssf}
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
/rlt/SetFileName root_file/mup_{muon_energy:.1f}GeV_Zp_{zp_mass:.2f}GeV_mumu_x{xssf:.0e}.root
/run/printProgress 1000
/run/beamOn 100000
'''[1:]

# Scan kinematics points.
ki_point_index = { }  # [zp_mass][muon_energy]
for ki_point in ki_points:
    ki_point_match = ki_point_pattern.search(ki_point)
    if not ki_point_match: continue
    zp_mass = float(ki_point_match.group(1))
    ki_point_index[zp_mass] = ki_point
ki_points = sorted(ki_point_index.items())

for zp_mass, ki_point in ki_points:

    # Compute cross section scale factor.
    with uproot.open(ki_point) as zp_mass_points:
        zp_mass_points = zp_mass_points['ki_points'].arrays(['zp_mass', 'muon_energy', 'xs'])
    assert np.all(zp_mass_points['zp_mass'] == zp_mass)
    spline = make_interp_spline(zp_mass_points['muon_energy'], zp_mass_points['xs'])
    xses = spline(muon_energies)
    xssfs = 10**np.ceil(np.log10(target_xs / xses))

    # Generate mac file.
    for muon_energy, xssf in zip(muon_energies, xssfs):
        mac_content = mac_template.format(ki_point='../' + ki_point, xssf=xssf, muon_energy=muon_energy, zp_mass=zp_mass)
        mac_path = 'build/mup_{muon_energy:.1f}GeV_Zp_{zp_mass:.2f}GeV_mumu_x{xssf:.0e}.mac'.format(
                muon_energy=muon_energy, zp_mass=zp_mass, xssf=xssf)
        with open(mac_path, 'w') as mac_file: mac_file.write(mac_content)
