"""
Step 0 validation for hybrid-path field extraction plan.

Compares the re-baselined 1×2 MMI hybrid S-matrix (current code, 2026-08-31)
against the known baseline from resolved/mmi_1x2_validation.md (2026-08-01).

Expected to match within numerical precision (both dispersive, both calibrated).
Any significant deviation indicates cache staleness or an unintended code change.
"""

import os
import pickle
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
duty_root = os.path.abspath(os.path.join(_HERE, "..", "..", "..", ".."))

# Load re-baselined results
rebaseline_pkl = os.path.join(duty_root, "example_output", "mmi_1x2", "mmi_1x2_hybrid", "smatrix_hybrid.pkl")
with open(rebaseline_pkl, 'rb') as f:
    smat_rebaseline = pickle.load(f)

# Extract transmission for both ports
wavelengths = smat_rebaseline.wavelengths
T_top = smat_rebaseline.transmission(port="out_top")
T_bot = smat_rebaseline.transmission(port="out_bot")

print("=" * 70)
print("STEP 0 COMPARISON: 1×2 MMI Hybrid Re-baseline")
print("=" * 70)
print(f"\nRe-baselined results from: {rebaseline_pkl}")
print(f"Wavelengths: {wavelengths[0]:.3f} – {wavelengths[-1]:.3f} µm ({len(wavelengths)} pts)")
print(f"Configuration: res=40, N_MODES=8, calibrate_body_neff=True, dispersive=True")
print()

# Baseline numbers from mmi_1x2_validation.md gap 5 (dispersive, calibrated)
# These are the FDTD reference values; the hybrid is expected to be ~10% lower
baseline_fdtd_t = np.array([
    0.44716,  # 1.260
    0.44068,  # 1.270
    0.42384,  # 1.280
    0.41259,  # 1.290
    0.40821,  # 1.300
    0.39527,  # 1.310 — center wavelength
    0.38649,  # 1.320
    0.37839,  # 1.330
    0.35202,  # 1.340 — large shift
    0.33743,  # 1.350
    0.34845,  # 1.360 — band edge (worst point)
])

print("Wavelength-by-wavelength comparison (hybrid transmission vs FDTD reference):")
print()
print(f"{'λ (µm)':<8} {'Hybrid T':<10} {'FDTD Ref':<10} {'ΔT':<10} {'Rel %':<10} {'Notes':<20}")
print("-" * 80)

delta_t = T_top - baseline_fdtd_t
rel_pct = 100.0 * delta_t / baseline_fdtd_t

for i, (wl, hy, fd, dt, rp) in enumerate(zip(wavelengths, T_top, baseline_fdtd_t, delta_t, rel_pct)):
    note = ""
    if i == len(wavelengths) // 2:
        note = "center"
    elif i == len(wavelengths) - 1:
        note = "band edge (worst)"
    elif wl == 1.34:
        note = "large shift"

    print(f"{wl:<8.3f} {hy:<10.5f} {fd:<10.5f} {dt:+10.5f} {rp:+10.2f}% {note:<20}")

print()
print("Summary statistics:")
print("-" * 80)
mean_abs_delta_t = np.mean(np.abs(delta_t))
max_abs_delta_t = np.max(np.abs(delta_t))
mean_t = np.mean(T_top)
rel_error = 100.0 * mean_abs_delta_t / np.mean(baseline_fdtd_t)

print(f"Mean |ΔT|:              {mean_abs_delta_t:.5f}")
print(f"Max  |ΔT|:              {max_abs_delta_t:.5f} at λ={wavelengths[np.argmax(np.abs(delta_t))]:.3f} µm")
print(f"Mean T (hybrid):        {mean_t:.5f}")
print(f"Mean T (FDTD):          {np.mean(baseline_fdtd_t):.5f}")
print(f"Relative error:         {rel_error:.2f}%")
print()

# Check against baseline expectation: 10.18% (mean rel error from 2026-08-01)
baseline_rel_error = 10.18
tolerance = 0.5  # allow ±0.5% deviation

print("Baseline expectation (from mmi_1x2_validation.md gap 5, 2026-08-01):")
print(f"  Mean |ΔT|:  0.02541")
print(f"  Rel error:  10.18%")
print()

if abs(rel_error - baseline_rel_error) < tolerance:
    print(f"✅ PASS: Relative error {rel_error:.2f}% matches baseline {baseline_rel_error:.2f}% (within {tolerance}%)")
    print("   Cache is clean and code is unchanged.")
else:
    print(f"⚠️  DEVIATION: Relative error {rel_error:.2f}% differs from baseline {baseline_rel_error:.2f}%")
    print(f"   Deviation: {rel_error - baseline_rel_error:+.2f}%")
    print()
    print("   Possible causes:")
    print("   1. Cache staleness (excite-normalization fix 2026-08-09 not auto-invalidated)")
    print("   2. Unintended code change since 2026-08-01")
    print("   3. Numerical variation (should be <0.1%)")
    print()
    print("   Next steps:")
    print("   - Check MPB seed log for near-degenerate excited modes")
    print("   - If found, clear cache and re-run")
    print("   - If no near-degenerate modes, investigate code changes")

print()
print("=" * 70)

# Also check splitter symmetry
T_bot_mean = np.mean(T_bot)
T_top_mean = np.mean(T_top)
asymmetry = 100.0 * abs(T_bot_mean - T_top_mean) / ((T_bot_mean + T_top_mean) / 2)
print(f"\nSplitter symmetry check:")
print(f"  Mean T (out_top): {T_top_mean:.5f}")
print(f"  Mean T (out_bot): {T_bot_mean:.5f}")
print(f"  Asymmetry: {asymmetry:.4f}%")
if asymmetry < 0.1:
    print(f"  ✅ PASS: Splitter is symmetric")
else:
    print(f"  ⚠️  WARNING: Splitter shows {asymmetry:.2f}% asymmetry (expected <0.1%)")
