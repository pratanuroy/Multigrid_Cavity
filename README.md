# Multigrid_Cavity
This program solves incompressible Navier Stokes equations for a 2D lid driven cavity flow using a non-linear multigrid method. 

The details of the multigrid method can be found in the following paper:

Roy, Pratanu, N. K. Anand, and Diego Donzis. "A Parallel Multigrid Finite-Volume Solver on a Collocated Grid for Incompressible Navier-Stokes Equations." Numerical Heat Transfer, Part B: Fundamentals 67.5 (2015): 376-409.

Link (Published article) : http://www.tandfonline.com/doi/abs/10.1080/10407790.2014.985980

Link (Accepted Manuscript): https://www.researchgate.net/profile/Pratanu_Roy/publication/272947323_A_Parallel_Multigrid_Finite-Volume_Solver_on_a_Collocated_Grid_for_Incompressible_Navier-Stokes_Equations/links/553050590cf27acb0de85ac5/A-Parallel-Multigrid-Finite-Volume-Solver-on-a-Collocated-Grid-for-Incompressible-Navier-Stokes-Equations.pdf

Notes on how to use the code:

The input parameters can be changed in the input_data.txt file. The grid size should be a multiple of 2, e.g. 128, 256, and so on. 
The number of multigrid levels (num_level) should be increased by 1 with successive increase in grid size. 

For postprocessing, set the write_vtk flag to .true. and the output will be written in VTR/VTK format. The .visit file can be loaded in Visit.
