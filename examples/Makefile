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

EXAMPLES = $(MMI_1X2_EXAMPLES) $(MMI_2X2_EXAMPLES) $(DIRECTIONAL_COUPLER_EXAMPLES) $(MODAL_MULTIPLEXER_EXAMPLES)
DEVICES = mmi_1x2 mmi_2x2 directional_coupler modal_multiplexer

# The brute-force FDTD coupling-length sweep is deliberately NOT in EXAMPLES
# or the directional_coupler aggregate: it's a long production sweep with its
# own batch script (run_directional_coupler_fdtd_sweep.batch, one SLURM array
# task per length), and keeping its name out of every other script's list is
# what guarantees no other batch script's per-example cleanup ever removes
# its output folder. SWEEP_INDEX selects one length (0..30); empty runs the
# whole sweep serially.
SWEEP_INDEX ?=

.PHONY: all $(DEVICES) $(EXAMPLES) directional_coupler_fdtd_sweep clean clean_directional_coupler_fdtd_sweep

all: $(EXAMPLES)

# Per-device aggregate targets (sequential; the per-device batch scripts run
# the same lists concurrently via SLURM job arrays instead).
mmi_1x2: $(MMI_1X2_EXAMPLES)
mmi_2x2: $(MMI_2X2_EXAMPLES)
directional_coupler: $(DIRECTIONAL_COUPLER_EXAMPLES)
modal_multiplexer: $(MODAL_MULTIPLEXER_EXAMPLES)

mmi_1x2_eme:
	$(runner) $(my_python) mmi_1x2/mmi_1x2_eme.py

mmi_1x2_fdtd:
	$(runner) $(my_python) mmi_1x2/mmi_1x2_fdtd.py

mmi_1x2_hybrid:
	$(runner) $(my_python) mmi_1x2/mmi_1x2_hybrid.py

mmi_2x2_eme:
	$(runner) $(my_python) mmi_2x2/mmi_2x2_eme.py

mmi_2x2_fdtd:
	$(runner) $(my_python) mmi_2x2/mmi_2x2_fdtd.py

mmi_2x2_hybrid:
	$(runner) $(my_python) mmi_2x2/mmi_2x2_hybrid.py

directional_coupler_hybrid:
	$(runner) $(my_python) directional_coupler/directional_coupler_hybrid.py

directional_coupler_mode:
	$(my_python) directional_coupler/directional_coupler_mode.py

directional_coupler_fdtd:
	$(runner) $(my_python) directional_coupler/directional_coupler_fdtd.py

directional_coupler_fdtd_sweep:
	$(runner) $(my_python) directional_coupler/directional_coupler_fdtd_sweep.py $(SWEEP_INDEX)

modal_multiplexer_hybrid:
	$(runner) $(my_python) modal_multiplexer/modal_multiplexer_hybrid.py

clean:
	rm -rf ../../example_output

clean_directional_coupler_fdtd_sweep:
	rm -rf ../../example_output/directional_coupler/directional_coupler_fdtd_sweep
