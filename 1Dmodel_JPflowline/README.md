Input files for running 1D idealised temperature simulations representing Lambert Amery glacial system

Notes:

To create a 1D mesh for simple temperature simulations: \
ElmerGrid 1 5 ColumnMesh.grd \
ElmerGrid 1 2 ColumnMesh.grd

Steady simulation is run by ColumnTemperature.sif.  Transient simulation (restarting
from steady simulation) is run be ColumnTemperature_transient.sif.
 
First line of ELMERSOLVER_STARTINFO determines which simulations is run (which .sif to use).

The .lua files contain some physical constants and functions.

Solver "saveline" outputs a text file (profile.dat) containing a vertical profile.
See profile.dat.names for variable names.
Example matlab script for plotting these is PlotProfiles.m.

Solver "ResultOutputSolver" writes the usual .vtu files.

See "Vertical Velocity" initial condition in steady .sif for setting vertical velocity
profile, with upper surface velocity set to approximate SMB.

See FS lower Accumulation Flux 3 for setting BMB in transient .sif.

