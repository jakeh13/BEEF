seed_num ?= 1
proc ?= 16
ifeq ($(LMOD_SYSHOST), phoenix)
	runner = srun -n
	trash = ~/scratch
	my_python = python3
else
	runner = mpirun -np
	trash = .
	my_python = python
endif

pol_demux_gc:
	$(runner) $(proc) $(my_python) src/pol_demux_grating.py -c my_configs/config_pol_demux_GC_$(seed_num).json -lv DEBUG -r opt

mmi:
	$(my_python) mmi_1x2_master.py -c configs_mmi/config.json

config:
	$(my_python) config.py
