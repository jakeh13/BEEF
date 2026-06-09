import os
import numpy as np
import meow as mw
from matplotlib import pyplot as plt

def plot_eme_fields(hp, cells, modes, tag=1, input_mode=1):
    fields_dir = f"{hp.outputfolder}/fields"
    os.makedirs(fields_dir, exist_ok=True)

    z = np.linspace(cells[0].z_min, cells[-1].z_max, 500)
    y = 0.0
    exc = input_mode - 1  # convert 1-based config index to 0-based meow index

    for component, label in [("Ex", "|Ex|"), ("Ey", "|Ey|"), ("Ez", "|Ez|"),
                              ("Hx", "|Hx|"), ("Hy", "|Hy|"), ("Hz", "|Hz|")]:
        field, x = mw.propagate_modes(modes, cells, y=y, z=z,
                                      excite_mode_l=exc, component=component)
        fig, ax = plt.subplots(figsize=(10, 4))
        im = ax.pcolormesh(z, x, np.abs(field).T, shading="auto", cmap="jet")
        ax.set_xlabel("z [um]")
        ax.set_ylabel("x [um]")
        ax.set_title(f"{label}(x, z) — input mode {input_mode}")
        plt.colorbar(im, ax=ax, label=label)
        plt.tight_layout()
        plt.savefig(f"{fields_dir}/eme_{component}_{tag}.png", dpi=300)
        plt.close("all")
