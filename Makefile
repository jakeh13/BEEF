seed_num ?= 1
ifeq ($(LMOD_SYSHOST), phoenix)
	runner = srun
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

config:
	$(my_python) config.py

config_cluster:
	$(runner) 1 $(my_python) config.py
