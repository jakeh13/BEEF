"""Determines which real "DC_test" GDS device id (1901/1902, 1911/1912,
1921/1922, 1931/1932 on the 4C10 die) corresponds to which coupling_length,
by comparing MEASURED spectral shape rather than re-deriving each one from
the die GDS directly (see directional_coupler_soi160_hybrid_dc_test.py's
own module docstring for why the GDS-geometry approach alone wasn't fully
reliable here -- 4 devices packed close together, plus heater/electrode
comb structures on the same GDS layer producing false rectangle-pair
candidates).

Data source (read-only; OneDrive, not part of this repo):
    .../Shared Documents - Ralph Group/Projects/Tapeouts/GF45CLO_4C10/
    4C10_measurements/{right_chip_0,right_chip_1}/*_DC_test_IL_sweep_*.json

Method: for each device pair, load both ids' AIN0/AIN1 channels (loopback-
normalized by id1012's own alignment_GC trace), and compute the mean
|AIN0-AIN1| dB split across the C-band. At a fixed gap (=0.25 um, the real
duty/dc.gds value), coupling ratio follows sin^2(kappa*L) -- the SMALLER
this split, the closer coupling_length sits to a coupling PEAK (near 50:50
through/cross); the LARGER the split, the closer to a coupling trough (weak
coupling). Devices are then ordered by this split and matched against the
3 already-known/hypothesized coupling_length candidates (0, 0.5, 11.634 um)
plus the one ALREADY simulated (3.0 um, directional_coupler_soi160_
hybrid.py) -- 1911/1912's own measured split (~10 dB) matches that already-
simulated device's own measured split (9.3-9.4 dB, from the LUNA-id
3913/3914/3923/3924 data in manuscript_figures/fp_resonator_figures/
HCG_plotting_LUNA.ipynb) closely enough that 1911/1912 is almost certainly
that SAME physical device under its original fab-time GDS label, not a
5th independent one.

Run:
    python /path/to/dc_test_id_to_coupling_length.py
"""

import json
from pathlib import Path

import numpy as np

DATA_DIR = Path(
    "/Users/jhiesener/Library/CloudStorage/OneDrive-GeorgiaInstituteofTechnology/"
    "Shared Documents - Ralph Group/Projects/Tapeouts/GF45CLO_4C10/4C10_measurements"
)
FOLDER = DATA_DIR / "right_chip_1"

DEVICE_PAIRS = {
    "1901/1902": (1901, 1902),
    "1911/1912": (1911, 1912),
    "1921/1922": (1921, 1922),
    "1931/1932": (1931, 1932),
}

# Coupling_length candidates -- 3.0 already simulated & matched to 1911/1912
# below; the other 3 are the working hypothesis this script's output backs.
CANDIDATE_LENGTHS_UM = [11.634, 3.0, 0.5, 0.0]


def load_id(folder: Path, dev_id: int) -> dict:
    matches = list(folder.glob(f"*id{dev_id}_*IL_sweep*.json"))
    if len(matches) != 1:
        raise FileNotFoundError(f"expected exactly 1 IL_sweep match for id{dev_id} in {folder}, got {matches}")
    with open(matches[0]) as f:
        return json.load(f)


def mean_split_db(folder: Path, id_a: int, id_b: int, lpbk_ain0: np.ndarray, wl: np.ndarray) -> float:
    """Mean |AIN0-AIN1| dB split across the C-band, averaged over both ids
    in the pair (each id excites a different input arm; by reciprocity the
    two should roughly agree -- averaging both is a cheap consistency
    check, not load-bearing for the ordering itself)."""
    mask = (wl >= 1525) & (wl <= 1600)
    splits = []
    for dev_id in (id_a, id_b):
        d = load_id(folder, dev_id)
        ain0 = np.array(d["values"]["transmission AIN0 [dBm]"]) - lpbk_ain0
        ain1 = np.array(d["values"]["transmission AIN1 [dBm]"]) - lpbk_ain0
        splits.append(np.abs(ain0[mask] - ain1[mask]).mean())
    return float(np.mean(splits))


def main():
    lpbk = load_id(FOLDER, 1012)
    wl = np.array(lpbk["values"]["wavelength [nm]"])
    lpbk_ain0 = np.array(lpbk["values"]["transmission AIN0 [dBm]"])

    results = {}
    for name, (id_a, id_b) in DEVICE_PAIRS.items():
        results[name] = mean_split_db(FOLDER, id_a, id_b, lpbk_ain0, wl)

    print("Device pair -> mean AIN0/AIN1 split (dB), C-band 1525-1600nm:")
    for name, split in sorted(results.items(), key=lambda kv: kv[1]):
        print(f"  {name:12s} {split:6.2f} dB")

    ordered = sorted(results.items(), key=lambda kv: kv[1])  # smallest split first
    print("\nSmallest split = closest to a coupling PEAK = LONGEST coupling_length"
          " (within this candidate set). Matched mapping:")
    for (name, split), length in zip(ordered, CANDIDATE_LENGTHS_UM):
        print(f"  {name:12s} split={split:6.2f} dB  ->  coupling_length = {length} um")


if __name__ == "__main__":
    main()
