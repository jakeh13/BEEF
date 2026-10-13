PYTHON ?= python
NP     ?= 4  # MPI processes for the non-Phoenix (mpirun) fallback only

ifeq ($(LMOD_SYSHOST), phoenix)
	runner = srun --mpi=pmi2
	my_python = python3
else
	runner = mpirun -np $(NP)
	my_python = $(PYTHON)
endif

# One folder per device under examples/, mirrored in ../../example_output/.
# Every device ships the same three methods where they are affordable:
# EME (fast, no FDTD), whole-device FDTD (the reference), and the hybrid
# cascade. Sweep drivers, verification harnesses, platform variants and
# retired examples moved to github.com/jakeh13/BEEF_archive on 2026-09-18;
# their targets went with them, which is why this file is much shorter than
# its history suggests.
MMI_1X2_EXAMPLES = mmi_1x2_eme mmi_1x2_fdtd mmi_1x2_hybrid
MMI_2X2_EXAMPLES = mmi_2x2_eme mmi_2x2_fdtd mmi_2x2_hybrid
# The WAISTED 2x2 (VALIDATED 2026-09-14). Its own device folder,
# examples/mmi_2x2_tapered/, holding the same three methods as the straight
# device. Target names carry the `_tapered_` infix because make targets share
# one namespace, while the scripts and their example_output/ subfolders are
# named for the METHOD (mmi_2x2_eme / _fdtd / _hybrid) inside the per-device
# folder, exactly as the straight device's are.
# ⚠ the HYBRID entry here is the SHARDED one, deliberately. Running that
# script once over all 11 wavelengths is a full cache MISS (~15 h): the
# solve_section cache key contains the whole wavelength array, so an 11-wl
# call and eleven 1-wl calls are different entries, and only the latter reuse
# what example_sim_cache/ already holds. `make all` must not silently cost 15
# hours, so the aggregate targets go through mmi_2x2_tapered_hybrid_all.
MMI_2X2_TAPERED_EXAMPLES = mmi_2x2_tapered_eme mmi_2x2_tapered_fdtd mmi_2x2_tapered_hybrid_all
MMI_1X4_EXAMPLES = mmi_1x4_eme mmi_1x4_fdtd mmi_1x4_hybrid
MMI_4X4_EXAMPLES = mmi_4x4_90hybrid_eme mmi_4x4_90hybrid_hybrid
DIRECTIONAL_COUPLER_EXAMPLES = directional_coupler_fdtd directional_coupler_hybrid
RING_RESONATOR_EXAMPLES = ring_resonator_fdtd ring_resonator_hybrid
HCG_EXAMPLES = high_contrast_grating_fdtd hcg_fabry_perot_hybrid

EXAMPLES = $(MMI_1X2_EXAMPLES) $(MMI_2X2_EXAMPLES) $(MMI_2X2_TAPERED_EXAMPLES) \
           $(MMI_1X4_EXAMPLES) $(MMI_4X4_EXAMPLES) $(DIRECTIONAL_COUPLER_EXAMPLES) \
           $(RING_RESONATOR_EXAMPLES) $(HCG_EXAMPLES)
DEVICES = mmi_1x2 mmi_2x2 mmi_2x2_tapered mmi_1x4 mmi_4x4_90hybrid \
          directional_coupler ring_resonator high_contrast_grating

.PHONY: all $(DEVICES) $(EXAMPLES) mmi_2x2_tapered_hybrid mmi_2x2_tapered_hybrid_merge clean

all: $(EXAMPLES)

# Per-device aggregate targets (sequential; the per-device batch scripts run
# the same lists concurrently via SLURM job arrays instead).
mmi_1x2: $(MMI_1X2_EXAMPLES)
mmi_2x2: $(MMI_2X2_EXAMPLES)
mmi_2x2_tapered: $(MMI_2X2_TAPERED_EXAMPLES)
mmi_1x4: $(MMI_1X4_EXAMPLES)
mmi_4x4_90hybrid: $(MMI_4X4_EXAMPLES)
directional_coupler: $(DIRECTIONAL_COUPLER_EXAMPLES)
ring_resonator: $(RING_RESONATOR_EXAMPLES)
high_contrast_grating: $(HCG_EXAMPLES)

# --- 1x2 MMI splitter -----------------------------------------------------
mmi_1x2_eme:
	$(runner) $(my_python) mmi_1x2/python/mmi_1x2_eme.py

mmi_1x2_fdtd:
	$(runner) $(my_python) mmi_1x2/python/mmi_1x2_fdtd.py

mmi_1x2_hybrid:
	$(runner) $(my_python) mmi_1x2/python/mmi_1x2_hybrid.py

# --- 2x2 MMI coupler ------------------------------------------------------
mmi_2x2_eme:
	$(runner) $(my_python) mmi_2x2/python/mmi_2x2_eme.py

mmi_2x2_fdtd:
	$(runner) $(my_python) mmi_2x2/python/mmi_2x2_fdtd.py

mmi_2x2_hybrid:
	$(runner) $(my_python) mmi_2x2/python/mmi_2x2_hybrid.py

# --- waisted 2x2 ----------------------------------------------------------
# Pure EME: no FDTD anywhere, so it is the cheap one -- minutes, against hours
# for the whole-device FDTD. Its body is NINE regions
# (wg|taper|stub|ramp|waist|ramp|stub|taper|wg) rather than the straight
# device's five, because a waisted body is three different things; the layout
# is the one the step-6 design sweep scored the grid on.
mmi_2x2_tapered_eme:
	$(runner) $(my_python) mmi_2x2_tapered/python/mmi_2x2_eme.py

mmi_2x2_tapered_fdtd:
	$(runner) $(my_python) mmi_2x2_tapered/python/mmi_2x2_fdtd.py

# The hybrid runs ONE WAVELENGTH PER INVOCATION and then merges: MPB does not
# scale with MPI ranks (meep solves the eigenproblem redundantly on every
# rank), so the parallelism has to live at the job level. MMI2X2T_HY_WLIDX
# selects the shard; MMI2X2T_HY_MERGE=1 combines them once all 11 have landed.
# See batch/run_mmi_2x2_tapered_examples.batch, which drives the whole thing
# as a SLURM array (task 0 = FDTD, tasks 1-11 = one wavelength each).
MMI2X2T_HY_WLIDX ?=
mmi_2x2_tapered_hybrid:
	MMI2X2T_HY_WLIDX=$(MMI2X2T_HY_WLIDX) $(runner) $(my_python) mmi_2x2_tapered/python/mmi_2x2_hybrid.py

mmi_2x2_tapered_hybrid_merge:
	MMI2X2T_HY_MERGE=1 $(my_python) mmi_2x2_tapered/python/mmi_2x2_hybrid.py

# Every wavelength, then the merge -- the whole hybrid leg in one target, and
# the one that `make mmi_2x2_tapered` / `make all` actually use. Sequential
# here because a plain `make` has no job scheduler; on the cluster, submit
# batch/run_mmi_2x2_tapered_examples.batch instead, which runs the same 11
# wavelengths as concurrent array tasks.
MMI2X2T_HY_NWL ?= 11
mmi_2x2_tapered_hybrid_all:
	@k=0; while [ $$k -lt $(MMI2X2T_HY_NWL) ]; do \
		$(MAKE) --no-print-directory mmi_2x2_tapered_hybrid MMI2X2T_HY_WLIDX=$$k || exit 1; \
		k=$$((k+1)); \
	done
	@$(MAKE) --no-print-directory mmi_2x2_tapered_hybrid_merge

# --- 1x4 MMI splitter -----------------------------------------------------
# mmi_1x4_fdtd was excluded here for a long time on the grounds that
# mm_length=206 um (216 with tapers) made a single-domain whole-device FDTD
# impractical. That reasoning went stale on 2026-08-14 when the body was
# re-optimised to mm_length=39.90 um, and the script was revived 2026-08-19.
mmi_1x4_eme:
	$(runner) $(my_python) mmi_1x4/python/mmi_1x4_eme.py

mmi_1x4_fdtd:
	$(runner) $(my_python) mmi_1x4/python/mmi_1x4_fdtd.py

mmi_1x4_hybrid:
	$(runner) $(my_python) mmi_1x4/python/mmi_1x4_hybrid.py

# --- 4x4 90-degree hybrid -------------------------------------------------
# These three had no make targets at all until 2026-09-18, which is the
# recorded reason the device went unrun for so long.
mmi_4x4_90hybrid_eme:
	$(runner) $(my_python) mmi_4x4_90hybrid/python/mmi_4x4_90hybrid_eme.py

mmi_4x4_90hybrid_hybrid:
	$(runner) $(my_python) mmi_4x4_90hybrid/python/mmi_4x4_90hybrid_hybrid.py

# --- directional coupler --------------------------------------------------
directional_coupler_fdtd:
	$(runner) $(my_python) directional_coupler/python/directional_coupler_fdtd.py

directional_coupler_hybrid:
	$(runner) $(my_python) directional_coupler/python/directional_coupler_hybrid.py

# --- ring resonator -------------------------------------------------------
# ⚠ the two ring scripts take DIFFERENT knobs. ring_resonator_fdtd.py does not
# read RING_TEST_MODE at all: it takes RING_FULL_3D (default 0, i.e. 2-D) and
# RING_FULL_RES (default 40). 2-D at res 40 is ~1.7 h; 3-D at the same
# resolution is ~167 h, which is why 2-D is the default and why this target
# does not force 3-D. Pass RING_FULL_3D=1 yourself if you have the allocation.
ring_resonator_fdtd:
	$(runner) $(my_python) ring_resonator/python/ring_resonator_fdtd.py

ring_resonator_hybrid:
	RING_TEST_MODE=0 $(runner) $(my_python) ring_resonator/python/ring_resonator_hybrid.py

# --- high-contrast grating ------------------------------------------------
# FDTD only (FDTDSolver.solve_reflection() — no EME, no HybridSolver).
# HCG_TEST_MODE=0 forces the full-res 3-D production run (defaults to the
# fast 2-D smoke test otherwise — see the script's own TEST_MODE docstring).
high_contrast_grating_fdtd:
	HCG_TEST_MODE=0 $(runner) $(my_python) high_contrast_grating/python/high_contrast_grating_fdtd.py

# The HCG mirror + directional coupler + connecting waveguides, assembled via
# beef.network.Circuit. It SOLVES all three pieces itself against
# example_sim_cache/ -- the mirror's solve_reflection(), the coupler's two
# solve_section() S-bends plus its beta calibration, and an FDTD beta
# calibration of the connecting waveguide. Cache hits make a re-run cheap, but
# a cold call is real 3-D FDTD and needs the ranks.
hcg_fabry_perot_hybrid:
	$(runner) $(my_python) high_contrast_grating_fabry_perot/python/hcg_fabry_perot_hybrid.py

clean:
	rm -rf ../../example_output
