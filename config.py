from engine.hyperparams.MMI import HyperParamsMMI, GeoParamsMMI, SimParams
from dataclasses import asdict
import json
import os
import numpy as np

OUTPUTFOLDER = "./out"
SCRATCH = "/storage/home/hcoda1/9/akhurana9/scratch/happi"

def save(hp: HyperParamsMMI, config_folder: str, name: str) -> None:
    s = json.dumps(asdict(hp), indent=4, default=lambda x: x.name)
    os.makedirs(config_folder, exist_ok=True)
    with open(f"{config_folder}/{name}.json", "w") as f:
        print(s, file=f)

def create_config_MMI_1x2(folder: str = 'configs_mmi', name: str = 'config'):
    res = 10
    sim_type = "EME"
    hp = HyperParamsMMI(
        geo=GeoParamsMMI(
            taper_width=0.8
        ),
        sim=SimParams(
            res=res,
            wavelength_min=1.27,
            wavelength_max=1.33,
            wavelength_pts=3,
            sim_type = sim_type
        ),
        outputfolder=f"./out/mmi_1x2/{sim_type}_{int(res)}_2"
    )
    save(hp, folder, name)

def create_config_MMI_2x2(folder: str = 'configs_mmi', name: str = 'config'):
    taper_spacing = 0.2 # microns between input and output tapers
    mm_width = 2.8
    mm_length = 28
    taper_width = mm_width/2 - taper_spacing/2
    waveguide_offset = taper_width + taper_spacing
    sim_type = "HYBRID"
    res = 40
    hp = HyperParamsMMI(
        geo=GeoParamsMMI(
            waveguide_width=0.45, # NOTE: Needs to be fixed for EME simulation waveguides (gdsfactory cross-section issue)
            mm_width=mm_width,
            mm_length=mm_length,
            taper_length=4,
            waveguide_offset=waveguide_offset,
            taper_width=taper_width,
            waveguide_height=0.22
        ),
        sim=SimParams(
            res=res,
            wavelength_min=1.530,
            wavelength_max=1.565,
            wavelength_pts=5,
            sim_type=sim_type,
            mm_monitor_position=0.25,
            input_modes=[1],
            output_modes=[1],
        ),
        outputfolder=f"./out/mmi_2x2_sweep/{sim_type}_{int(res)}_1"
    )
    save(hp, folder, name)

def create_config_MMI_2x2_gds(folder: str = 'configs_mmi', name: str = 'config_2x2_gds'):
    taper_spacing = 0.2
    mm_width = 1.6
    mm_length = 10.8
    taper_width = mm_width/2 - taper_spacing/2
    waveguide_offset = taper_width + taper_spacing
    sim_type = "FDTD_GDS"
    res = 40
    splitter_name = "Splitter_60"
    hp = HyperParamsMMI(
        geo=GeoParamsMMI(
            waveguide_width=0.5,
            mm_width=mm_width,
            mm_length=mm_length,
            taper_length=4,
            waveguide_offset=waveguide_offset,
            taper_width=taper_width,
            waveguide_height=0.215,
            gds_file=f"gds/{splitter_name}.gds"
        ),
        sim=SimParams(
            res=res,
            wavelength_min=1.520,
            wavelength_max=1.580,
            wavelength_pts=30,
            sim_type=sim_type,
            mm_monitor_position=0.25,
            input_modes=[1],
            output_modes=[1],
        ),
        outputfolder=f"./out/mmi_2x2_gds/{sim_type}_{int(res)}_{splitter_name}"
    )
    save(hp, folder, name)

def create_config_MMI_2x2_pol_splitter(folder: str = 'configs_mmi', name: str = 'config_2x2_gds'):
    taper_spacing = 0.2
    mm_width = 2.1
    mm_length = 37.9
    taper_width = mm_width/2 - taper_spacing/2
    waveguide_offset = taper_width + taper_spacing
    sim_type = "EME"
    res = 40
    # splitter_name = "Splitter_60"
    hp = HyperParamsMMI(
        geo=GeoParamsMMI(
            waveguide_width=0.35,
            mm_width=mm_width,
            mm_length=mm_length,
            taper_length=2,
            waveguide_offset=waveguide_offset,
            taper_width=taper_width,
            waveguide_height=0.16
        ),
        sim=SimParams(
            res=res,
            wavelength_min=1.270,
            wavelength_max=1.330,
            wavelength_pts=5,
            input_modes=[1,2],
            output_modes=[1,2],
            sim_type=sim_type,
            mm_monitor_position=0.5,
            num_modes_MMI=50, # 20
            num_modes_taper=20
        ),
        outputfolder=f"./out/mmi_2x2_pol_splitter/v81_{sim_type}_{int(res)}"
    )
    save(hp, folder, name)

def create_config_MMI_1x4(folder: str = 'configs_mmi', name: str = 'config_1x4'):
    waveguide_offset = 1.0  # pitch between adjacent output waveguide centers
    mm_width = 4.0          # needs to span 4 outputs: ~4 * waveguide_offset
    mm_length = 7.3
    taper_width = 0.8
    sim_type = "EME"
    res = 40
    hp = HyperParamsMMI(
        geo=GeoParamsMMI(
            mm_width=mm_width,
            mm_length=mm_length,
            taper_width=taper_width,
            waveguide_offset=waveguide_offset
        ),
        sim=SimParams(
            res=res,
            wavelength_min=1.26,
            wavelength_max=1.36,
            wavelength_pts=5,
            sim_type=sim_type
        ),
        outputfolder=f"./out/mmi_1x4_sweep/{sim_type}_{int(res)}_5"
    )
    save(hp, folder, name)

if __name__ == '__main__':
    create_config_MMI_1x2(name='config')
    create_config_MMI_2x2(name='config_2x2')
    create_config_MMI_2x2_gds(name='config_2x2_gds')
    create_config_MMI_2x2_pol_splitter(name='config_2x2_pol_splitter')
    create_config_MMI_1x4(name='config_1x4')