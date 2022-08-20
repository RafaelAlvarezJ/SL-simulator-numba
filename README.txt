Streamlines simulator for waterflooding in a quadrant of five-spot pattern, 
considering two-dimension incompressible flow and no capillary pressure. Transport 
calculations are just IMPES calculations (read "Streamlines Simulation: Theory and 
Practice"), which are faster than conventional finite-difference simulators.

It uses a fixed dt and a natural Time of Flight discretization (given by Pollock´s 
algorithm). It is optimized with Numba´s @njit decorator, powered by several 
compilation options.

The simulator consists of 3 files, and an optional one (Units): 
- SL_Main.py
- SL_Auxiliares_numba.py (Auxliar functions for the main program)
- SL_Datos.py (Simulation Data)
- Units.py (library with unit convertions, in order to work with consistent units)

All the simulation data is written in SL_Datos.py, and the velocity field is read at
SL_Auxiliares_numba.py's beginning (if you want to use another velocity field, here 
it requires to write its file name in order to switch velocity fields).

The velocity fields are from an IMPES simulator. The lines are traced from the 
injector well and this can´t be changed yet, such as the injector well cell. 
Nevertheless, the producer can be placed in any cell. Also, the velocity field 
must be from an squared mesh (30x30, 60x60, 200x200, and so on).
