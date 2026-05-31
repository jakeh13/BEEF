from engine.hyperparams.MMI_1x2 import HyperParamsMMI, GeoParamsMMI, SimParams
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
    hp = HyperParamsMMI(
        geo=GeoParamsMMI(
            taper_width=0.8
        ),
        sim=SimParams(
            res=40,
            wavelength_min=1.27,
            wavelength_max=1.33,
            wavelength_pts=20,
            sim_type = "EME"
        ),
        outputfolder="./out/mmi_1x2_eme_1"
    )
    save(hp, folder, name)

def create_config_MMI_2x2(folder: str = 'configs_mmi', name: str = 'config'):
    taper_spacing = 0.2 # microns between input and output tapers
    mm_width = 2.8
    mm_length = 28
    taper_width = mm_width/2 - taper_spacing/2
    waveguide_offset = taper_width + taper_spacing
    sim_type = "HYBRID"
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
            res=20,
            wavelength_min=1.535,
            wavelength_max=1.565,
            wavelength_pts=5,
            sim_type = sim_type,
            mm_monitor_position=0.25
        ),
        outputfolder=f"./out/mmi_2x2_220/{sim_type}_{int(mm_length)}_1"
    )
    save(hp, folder, name)

if __name__ == '__main__':
    create_config_MMI_1x2(name='config')
    create_config_MMI_2x2(name='config_2x2')