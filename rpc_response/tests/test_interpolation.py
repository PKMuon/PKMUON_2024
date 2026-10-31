"""Behavior checks on real library tables, including validity propagation."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from response_interpolation import MuonStandardInterpolator

library = MuonStandardInterpolator(sys.argv[1])
a, b, w = library.bracket(np.sqrt(.3))
assert (a, b) == (.3, 1.) and abs(w-.5) < 1e-12
frame, weights = library.distribution(.6)
assert abs(weights.sum()-1) < 1e-12 and (weights >= 0).all()
sample = library.sample(.6, 10000, 123)
assert len(sample) == 10000
assert (sample.n_electrons % 1 == 0).all()
for source, group in sample.groupby("source_standard_energy_gev"):
    reference = library.frames[source].set_index("event").loc[group.source_standard_event]
    for field in ("n_clusters", "n_electrons", "energy_transfer_ev", "charge_valid", "charge_fc"):
        assert np.allclose(group[field], reference[field], equal_nan=True)
assert sample.loc[sample.charge_valid == 0, "charge_fc"].isna().all()
for energy in (0., .09, 101., float("nan")):
    try:
        library.sample(energy, 1)
    except ValueError:
        pass
    else:
        raise AssertionError("Out-of-domain query was not rejected")
held = MuonStandardInterpolator(sys.argv[1], exclude_gev=1.)
assert 1. not in held.frames and held.bracket(1.)[:2] == (.3, 3.)
print("PASS exact weights, energy bounds, held-out exclusion, event correlation and invalid-charge propagation")
