import numpy as np
import pickle
from matplotlib import pyplot as plt

OUTPUT_FOLDER = "/storage/scratch1/8/jhiesener3/shit/out/mmi_1x2_gds/FDTD_GDS_40_splitter_f32"

with open(f"{OUTPUT_FOLDER}/data/sim_results.pkl", "rb") as f:
    data = pickle.load(f)

wavelengths = np.array(data["wavelengths"])
powers_top = data["powers_top_by_mode"]
powers_bot = data["powers_bot_by_mode"]
keys = sorted(powers_top.keys())

fig, ax = plt.subplots()
for in_m, out_m in keys:
    total = np.array(powers_top[(in_m, out_m)]) + np.array(powers_bot[(in_m, out_m)])
    ax.plot(wavelengths, total, label=f"In {in_m} → Out {out_m}")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Total Transmission")
ax.legend()
fig.tight_layout()
fig.savefig(f"{OUTPUT_FOLDER}/total_transmission.png", dpi=300)
plt.close(fig)

fig, ax = plt.subplots()
for in_m, out_m in keys:
    total = np.array(powers_top[(in_m, out_m)]) + np.array(powers_bot[(in_m, out_m)])
    ax.plot(wavelengths, 10 * np.log10(total), label=f"In {in_m} → Out {out_m}")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Total Transmission (dB)")
ax.legend()
fig.tight_layout()
fig.savefig(f"{OUTPUT_FOLDER}/total_transmission_dB.png", dpi=300)
plt.close(fig)

fig, ax = plt.subplots()
for in_m, out_m in keys:
    top = np.array(powers_top[(in_m, out_m)])
    total = top + np.array(powers_bot[(in_m, out_m)])
    ax.plot(wavelengths, top / total, label=f"In {in_m} → Out {out_m}")
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Splitting Ratio (Top / Total)")
ax.legend()
fig.tight_layout()
fig.savefig(f"{OUTPUT_FOLDER}/splitting_ratio.png", dpi=300)
plt.close(fig)
