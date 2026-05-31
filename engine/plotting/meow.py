import numpy as np
import meow as mw
from matplotlib import pyplot as plt

def plot_eme_fields(hp,cells,modes):
    # Propagate the modes through the structure and visualize the field distribution
    z = np.linspace(cells[0].z_min, cells[-1].z_max, 500)
    y = 0.0
    field, x = mw.propagate_modes(modes, cells, y=y, z=z)
    fig, ax = plt.subplots(figsize=(10, 4))
    im = ax.pcolormesh(z, x, np.abs(field).T, shading="auto", cmap="jet")
    ax.set_xlabel("z [um]")
    ax.set_ylabel("x [um]")
    ax.set_title("|Ex(x, z)| from mw.propagate_modes")
    plt.colorbar(im, ax=ax, label="|Ex|")
    plt.tight_layout()
    plt.savefig(f'{hp.outputfolder}/eme_fields.png', dpi=300)
