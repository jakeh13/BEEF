seed_num ?= 1
DILATION_NM ?= 0
MM_LENGTH_X10 ?= 40
export PYTHONPATH := $(CURDIR)
ifeq ($(LMOD_SYSHOST), phoenix)
	runner = srun --mpi=pmi2
	trash = ~/scratch
	my_python = python3
else
	runner = mpirun -np 2 
	trash = .
	my_python = python
endif

mmi_cluster:
	$(runner) $(my_python) masters/mmi_1x2_master.py -c configs_mmi/config.json

mmi_2x2_bend:
	$(runner) $(my_python) masters/mmi_2x2_master.py -c configs_mmi/config_2x2_gds_bend.json

mmi:
	$(runner) $(my_python) masters/mmi_1x2_master.py -c configs_mmi/config.json

mmi_1x2_gds:
	$(runner) $(my_python) masters/mmi_1x2_master.py -c configs_mmi/config_1x2_gds.json

mmi_2x2:
	$(runner) $(my_python) masters/mmi_2x2_master.py -c configs_mmi/config_2x2.json

mmi_2x2_gds:
	$(runner) $(my_python) masters/mmi_2x2_master.py -c configs_mmi/config_2x2_gds.json

mmi_2x2_pol_splitter:
	$(runner) $(my_python) masters/mmi_2x2_master.py -c configs_mmi/config_2x2_pol_splitter.json

mmi_2x2_s_cluster::
	$(runner) $(my_python) sweeps/mmi_2x2_sweep.py -c configs_mmi/config_2x2.json

mmi_2x2_s_eme_cluster:
	$(runner) $(my_python) sweeps/mmi_2x2_sweep_efficient.py -c configs_mmi/config_2x2.json

mmi_1x4:
	$(runner) $(my_python) masters/mmi_1x4_master.py -c configs_mmi/config_1x4.json

mmi_1x4_s:
	$(runner) $(my_python) sweeps/mmi_1x4_sweep_efficient.py --config configs_mmi/config_1x4.json

mmi_1x4_s_simple:
	$(runner) $(my_python) sweeps/mmi_1x4_sweep.py -c configs_mmi/config_1x4.json

mmi_mux:
	$(runner) $(my_python) masters/mmi_mux_master.py -c configs_mmi/config_mux.json

mmi_2x2_dilation_sweep:
	$(runner) $(my_python) masters/mmi_2x2_master.py -c configs_mmi/config_2x2_gds_dilation_$(DILATION_NM)nm.json

config_dilation_sweep_cluster:
	$(runner) $(my_python) config.py --dilation $(DILATION_NM)

mmi_1x2_length_sweep:
	$(runner) $(my_python) masters/mmi_1x2_master.py -c configs_mmi/config_1x2_gds_length_$(MM_LENGTH_X10).json

mmi_1x2_s_simple:
	$(runner) $(my_python) sweeps/mmi_1x2_sweep.py -c configs_mmi/config_1x2_sweep.json 

mmi_1x2_s_parallel:
	$(my_python) sweeps/mmi_1x2_sweep_parallel.py -c configs_mmi/config_1x2_sweep.json

mmi_1x2_s_parallel_efficient:
	$(my_python) sweeps/mmi_1x2_sweep_parallel_efficient.py -c configs_mmi/config_1x2_sweep.json --full-output

mmi_2x2_s_parallel_efficient:
	$(my_python) sweeps/mmi_2x2_sweep_parallel_efficient.py -c configs_mmi/config_2x2_sweep.json --full-output

config_1x2_length_sweep_cluster:
	$(runner) $(my_python) config.py --length_x10 $(MM_LENGTH_X10)

config:
	$(my_python) config.py

config_cluster:
	$(runner) $(my_python) config.py
