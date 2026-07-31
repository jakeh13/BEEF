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
# Keep each device's list in sync with its run_<device>_examples.batch
# EXAMPLES array at the duty root (index order matters there since
# SLURM_ARRAY_TASK_ID selects a name positionally).
MMI_1X2_EXAMPLES = mmi_1x2_eme mmi_1x2_fdtd mmi_1x2_hybrid
MMI_2X2_EXAMPLES = mmi_2x2_eme mmi_2x2_fdtd mmi_2x2_hybrid
DIRECTIONAL_COUPLER_EXAMPLES = directional_coupler_hybrid directional_coupler_mode directional_coupler_fdtd
MODAL_MULTIPLEXER_EXAMPLES = modal_multiplexer_hybrid
# modal_multiplexer_fdtd (gap 3's full-device reference, added 2026-07-31) is
# deliberately NOT in this list: run_modal_multiplexer_examples.batch is a
# single-task `--array=0-0` positionally bound to index 0 = hybrid, and a job
# is running off that exact file as this comment is written. Given its own
# batch script (run_modal_multiplexer_fdtd.batch) instead of reshaping a
# live array.
RING_RESONATOR_EXAMPLES = ring_resonator_hybrid ring_resonator_hybrid_cascaded

EXAMPLES = $(MMI_1X2_EXAMPLES) $(MMI_2X2_EXAMPLES) $(DIRECTIONAL_COUPLER_EXAMPLES) $(MODAL_MULTIPLEXER_EXAMPLES) $(RING_RESONATOR_EXAMPLES)
DEVICES = mmi_1x2 mmi_2x2 directional_coupler modal_multiplexer ring_resonator

# The five directional-coupler sweep drivers (brute-force FDTD coupling-
# length, N_MODES, mesh-resolution x2, body-stub) moved to
# examples/directional_coupler/archive/ on 2026-07-31, retired pending new,
# more general sweeping functionality once every device is verified (see
# notes.md Misc #5). Their targets, SWEEP_INDEX, and clean_* rules were
# removed from here with them; batch/archive/ holds their old batch scripts.
# test_eme_convergence.py still reads the archived stub/modes sweeps directly
# by path to guard the item-11 resolution-knob regression they're each
# checked for, so they're not orphaned, just not runnable via `make` until
# revived.

# TEMPORARY — verification harness for dc_validation.md items 1-2, kept out of
# EXAMPLES for the same reason the sweeps are (nothing else's cleanup should
# touch its output). Delete this target, clean_verification, and
# examples/verification/ once those items are signed off.
VERIFICATION = verify_items_1_2 verify_cache_multirank verify_item_5_ground_truth verify_item_6_engine_agreement

.PHONY: all $(DEVICES) $(EXAMPLES) $(VERIFICATION) clean clean_verification

all: $(EXAMPLES)

# Per-device aggregate targets (sequential; the per-device batch scripts run
# the same lists concurrently via SLURM job arrays instead).
mmi_1x2: $(MMI_1X2_EXAMPLES)
mmi_2x2: $(MMI_2X2_EXAMPLES)
directional_coupler: $(DIRECTIONAL_COUPLER_EXAMPLES)
modal_multiplexer: $(MODAL_MULTIPLEXER_EXAMPLES)
ring_resonator: $(RING_RESONATOR_EXAMPLES)

mmi_1x2_eme:
	$(runner) $(my_python) mmi_1x2/python/mmi_1x2_eme.py

mmi_1x2_fdtd:
	$(runner) $(my_python) mmi_1x2/python/mmi_1x2_fdtd.py

mmi_1x2_hybrid:
	$(runner) $(my_python) mmi_1x2/python/mmi_1x2_hybrid.py

mmi_2x2_eme:
	$(runner) $(my_python) mmi_2x2/python/mmi_2x2_eme.py

mmi_2x2_fdtd:
	$(runner) $(my_python) mmi_2x2/python/mmi_2x2_fdtd.py

mmi_2x2_hybrid:
	$(runner) $(my_python) mmi_2x2/python/mmi_2x2_hybrid.py

directional_coupler_hybrid:
	$(runner) $(my_python) directional_coupler/python/directional_coupler_hybrid.py

directional_coupler_mode:
	$(my_python) directional_coupler/python/directional_coupler_mode.py

directional_coupler_fdtd:
	$(runner) $(my_python) directional_coupler/python/directional_coupler_fdtd.py

# TEMPORARY (see VERIFICATION above). VERIFY_ARGS=--static-only runs just the
# no-meep checks, which is safe on a login node.
VERIFY_ARGS ?=
verify_items_1_2:
	$(runner) $(my_python) verification/python/verify_items_1_2.py $(VERIFY_ARGS)

# Cache write-gate check (item 4). Needs >=2 NODES to mean anything -- a
# single-node run cannot race across hosts -- and save_dir must be on shared
# storage, so run it through srun directly rather than the $(runner) default:
#   srun -A <account> -q inferno -N 2 --ntasks-per-node=24 -t 10 \
#        --mem-per-cpu=2G --mpi=pmi2 make verify_cache_multirank
verify_cache_multirank:
	$(runner) $(my_python) verification/python/verify_cache_multirank.py \
		../../example_output/verification/cache_multirank

# Ground-truth FDTD length sweep (item 5). GT_INDEX selects one of the 40
# lengths (0..39), mirroring SWEEP_INDEX above; empty runs the whole sweep
# serially. VERIFY_ARGS=--aggregate-only re-aggregates without simulating,
# which is safe on a login node; VERIFY_ARGS=--no-cache forces re-simulation.
GT_INDEX ?=
verify_item_5_ground_truth:
	$(runner) $(my_python) verification/python/verify_item_5_ground_truth.py \
		$(GT_INDEX) $(VERIFY_ARGS)

# Direct FDTD/EME engine-agreement measurement (item 6). EA_INDEX selects one of
# the 4 section lengths (0..3); empty runs all four serially. The EME half runs
# inside the aggregate step and needs no MPI, so
# VERIFY_ARGS=--aggregate-only gives the full comparison from whatever FDTD
# lengths exist -- safe on a login node.
EA_INDEX ?=
verify_item_6_engine_agreement:
	$(runner) $(my_python) verification/python/verify_item_6_engine_agreement.py \
		$(EA_INDEX) $(VERIFY_ARGS)

modal_multiplexer_hybrid:
	$(runner) $(my_python) modal_multiplexer/python/modal_multiplexer_hybrid.py

modal_multiplexer_fdtd:
	$(runner) $(my_python) modal_multiplexer/python/modal_multiplexer_fdtd.py

# RING_TEST_MODE=0 forces the full-res 3-D production run (defaults to the
# fast 2-D smoke test otherwise — see the script's own TEST_MODE docstring).
ring_resonator_hybrid:
	RING_TEST_MODE=0 $(runner) $(my_python) ring_resonator/python/ring_resonator_hybrid.py

ring_resonator_hybrid_cascaded:
	RING_TEST_MODE=0 $(runner) $(my_python) ring_resonator/python/ring_resonator_hybrid_cascaded.py

clean:
	rm -rf ../../example_output

clean_verification:
	rm -rf ../../example_output/verification
