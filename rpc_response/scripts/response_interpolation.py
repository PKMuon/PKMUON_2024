"""Event-preserving stochastic interpolation of compatible muon standards.

P_E=(1-w) P_low+w P_high, w=log(E/E_low)/log(E_high/E_low).
Each draw selects an entire existing event, preserving integer multiplicities,
within-event correlations and invalid-charge flags. This is a distribution
approximation for normal-incidence muons, not an avalanche saturation model.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

CONDITIONS = ("schema", "angle_deg", "sigma_input_mm", "gas_gap_mm", "temperature_k",
              "pressure_torr", "gas_voltage_v", "velocity_mm_ns", "alpha_per_mm",
              "eta_per_mm", "heed_w_ev")


class MuonStandardInterpolator:
    def __init__(self, directory, exclude_gev=None):
        self.frames, self.metadata = {}, {}
        reference = None
        for path in sorted(Path(directory).glob("*_metadata.json")):
            meta = json.loads(path.read_text())
            energy = meta["muon_kinetic_ev"] / 1e9
            if exclude_gev is not None and np.isclose(energy, exclude_gev, rtol=1e-12):
                continue
            if reference is None:
                reference = meta
            if any(meta.get(k) != reference.get(k) for k in CONDITIONS):
                raise ValueError("Incompatible standard gas/geometry/readout conditions")
            frame = pd.read_csv(str(path).replace("_metadata.json", "_events.csv"))
            if len(frame) != meta["events"] or energy in self.frames:
                raise ValueError("Incomplete/duplicate standard")
            self.frames[energy], self.metadata[energy] = frame, meta
        self.energies = np.array(sorted(self.frames))
        if not len(self.energies):
            raise ValueError("No standards")

    def bracket(self, energy_gev, method="log_mixture"):
        e = float(energy_gev)
        if not np.isfinite(e) or e < self.energies[0] or e > self.energies[-1]:
            raise ValueError("Energy outside standard range; extrapolation is disabled")
        if method not in ("nearest", "log_mixture", "linear_mixture"):
            raise ValueError("Unknown interpolation method")
        if method == "nearest":
            q = float(self.energies[np.argmin(np.abs(np.log(self.energies / e)))])
            return q, q, 0.
        exact = self.energies[self.energies == e]
        if len(exact):
            return e, e, 0.
        index = np.searchsorted(self.energies, e)
        low, high = map(float, self.energies[index-1:index+1])
        w = (np.log(e/low)/np.log(high/low) if method == "log_mixture"
             else (e-low)/(high-low))
        return low, high, float(w)

    def distribution(self, energy_gev, method="log_mixture"):
        low, high, w = self.bracket(energy_gev, method)
        if low == high:
            frame = self.frames[low]
            return frame, np.full(len(frame), 1/len(frame))
        a, b = self.frames[low], self.frames[high]
        return pd.concat([a, b], ignore_index=True), np.r_[
            np.full(len(a), (1-w)/len(a)), np.full(len(b), w/len(b))]

    def sample(self, energy_gev, count, seed=20261021, method="log_mixture"):
        low, high, w = self.bracket(energy_gev, method)
        if not isinstance(count, (int, np.integer)) or count <= 0:
            raise ValueError("Sample count must be a positive integer")
        rng = np.random.default_rng(seed)
        sources = np.where(rng.random(count) < w, high, low)
        rows = np.empty(count, dtype=int)
        output = []
        for source in np.unique(sources):
            selected = np.flatnonzero(sources == source)
            rows[selected] = rng.integers(0, len(self.frames[source]), len(selected))
            part = self.frames[source].iloc[rows[selected]].copy()
            part["sample_index"] = selected
            part["source_standard_energy_gev"] = source
            part["source_standard_event"] = part["event"]
            output.append(part)
        result = pd.concat(output).sort_values("sample_index").reset_index(drop=True)
        result["query_energy_gev"] = energy_gev
        result["interpolation_method"] = method
        return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("standards"); parser.add_argument("energy_gev", type=float)
    parser.add_argument("output"); parser.add_argument("--count", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261021)
    args = parser.parse_args()
    path = Path(args.output)
    if path.exists():
        raise FileExistsError("Output protected")
    MuonStandardInterpolator(args.standards).sample(args.energy_gev, args.count, args.seed).to_csv(path, index=False)
