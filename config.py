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
            res=30,
            wavelength_min=1.26,
            wavelength_max=1.36,
            wavelength_pts=10,
            sim_type = "HYBRID"
        ),
        outputfolder="./out/test_hybrid"
    )
    save(hp, folder, name)

if __name__ == '__main__':
    create_config_MMI_1x2(name='config')