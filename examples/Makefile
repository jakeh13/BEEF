PYTHON ?= python
NP     ?= 4  # number of MPI processes: make NP=4 mmi_1x2_fdtd

RUN = $(if $(filter 1,$(NP)),$(PYTHON),mpirun -np $(NP) $(PYTHON))

.PHONY: all mmi_1x2_eme mmi_1x2_fdtd mmi_1x2_hybrid clean

all: mmi_1x2_eme mmi_1x2_fdtd mmi_1x2_hybrid

mmi_1x2_eme:
	$(RUN) mmi_1x2_eme.py

mmi_1x2_fdtd:
	$(RUN) mmi_1x2_fdtd.py

mmi_1x2_hybrid:
	$(RUN) mmi_1x2_hybrid.py

clean:
	rm -rf ../../example_output
