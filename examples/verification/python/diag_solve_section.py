import numpy as np
from beef import MMI, EMESolver
from beef.materials.stack import SIN_300

mmi = MMI(
    n_inputs=1, n_outputs=4,
    mm_length=150.0, mm_width=12.0,
    waveguide_width=1.1, taper_length=5.0, taper_width=2.2,
    waveguide_offset=2.4, waveguide_height=0.30,
)
solver = EMESolver(
    num_modes=[6], num_cells=[1], lengths=[1.0],
    oxide_width=20.0, oxide_height=5.0, res=20.0,
    material_stack=SIN_300, mode_solver="tidy3d",
)
half_body = 150.0/2
half_total = mmi.total_length/2
x0, x1 = half_body, half_total  # output taper+access band, same span used by M_out
wl = np.array([1.55])
M = solver.solve_section(mmi, wl, x0, x1, n_modes_in=6, n_modes_out=6, output_idx=0)
print("M_out-band matrix, |.|, wl=1.55:")
np.set_printoptions(precision=3, suppress=True)
print(np.abs(M[0]))
