seed_num ?= 1
ifeq ($(LMOD_SYSHOST), phoenix)
	runner = srun --mpi=pmi2
	trash = ~/scratch
	my_python = python3
else
	runner = mpirun -np
	trash = .
	my_python = python
endif

mmi_cluster:
	$(runner) $(my_python) mmi_1x2_master.py -c configs_mmi/config.json

mmi:
	$(runner) 1 $(my_python) mmi_1x2_master.py -c configs_mmi/config.json

mmi_2x2:
	$(runner) 4 $(my_python) mmi_2x2_master.py -c configs_mmi/config_2x2.json

mmi_2x2_s:
	$(runner) 4 $(my_python) mmi_2x2_sweep.py -c configs_mmi/config_2x2.json

mmi_2x2_s_eme:
	$(runner) 4 $(my_python) mmi_2x2_sweep.py -c configs_mmi/config_2x2.json

mmi_1x4:
	$(runner) 2 $(my_python) mmi_1x4_master.py -c configs_mmi/config_1x4.json

mmi_1x4_s:
	$(runner) 2 $(my_python) mmi_1x4_sweep_efficient.py --config configs_mmi/config_1x4.json

config:
	$(my_python) config.py

config_cluster:
	$(runner) 1 $(my_python) config.py
