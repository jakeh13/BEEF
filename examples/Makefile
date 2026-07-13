PYTHON ?= python
NP     ?= 4  # MPI processes for the non-Phoenix (mpirun) fallback only

ifeq ($(LMOD_SYSHOST), phoenix)
	runner = srun --mpi=pmi2
	my_python = python3
else
	runner = mpirun -np $(NP)
	my_python = $(PYTHON)
endif

.PHONY: all mmi_1x2_eme mmi_1x2_fdtd mmi_1x2_hybrid clean

all: mmi_1x2_eme mmi_1x2_fdtd mmi_1x2_hybrid

mmi_1x2_eme:
	$(runner) $(my_python) mmi_1x2_eme.py

mmi_1x2_fdtd:
	$(runner) $(my_python) mmi_1x2_fdtd.py

mmi_1x2_hybrid:
	$(runner) $(my_python) mmi_1x2_hybrid.py

clean:
	rm -rf ../../example_output
