FDTD/EME combined electromagnetic solver enables the use of the MEEP FDTD solver and MEOW EME solver in simulating MMI devices.

To-do:
 - Parallelize the mode calculation using MPI
 - Develop reduced mode computation method that limits the mode calculation to ones with moderate power (i.e. specify the computed modes)
 - Investigate method to reduce 2x2 MMI splitter to one FDTD simulation
 - Add directional coupler simulation capability
 - Add edge-coupling spot-size converter simulation capability
 - Add bragg grating simulation capability
