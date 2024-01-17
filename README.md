![im](https://github.com/FChauvet/fractal_growth_DBM/assets/123599610/d81bb8d3-309f-43a8-b7bf-d66af2082e0f)

fractal_growth_vg is a code which simulates the growth, limited by diffusion and without stabilizing effects, of a 2D fractal deposit from a Laplacian growth type model (also known as a dielectric break model) and taking into account a positive growth rate (vg).

The concentration field around the deposit is computed by taking into account the "convection" effect induced by the front motion. To do this, the code solves an advection diffusion equation, using an iterative method, with a finite difference scheme. Starting from a flat deposit surface, at each iteration, the algorithm randomly selects a site on the deposit surface (growth on a lattice) from the distribution of the local mass flux (from the calculation of the concentration field). The calculation stops when the maximum height of the deposit reaches the desired value (hmax).
The domain must be large enough so that growth is not influenced by the boundary conditions (lateral = periodic, upper edge = Dirichlet (~1), deposit surface = Dirichlet (0)). The upper edge is automatically positioned at a sufficiently high distance from the highest point of the deposit surface (domain_height parameter).
The concentration matrix is solved by a fortran script loaded by python using numpy.f2py. The produced data (concentration and deposit matrices, parameters, etc.) are saved in a hdf5 file format.

Requirements :
- python 3 with numpy, scipy, matplotlib and h5py
- linux or os with an intel proc (for f2py)


Files :
- fractal_growth_vg.py (main program)
- solve_adv_diff.f90 (fortran script for solving the advection diffusion equation)
- script_showlastim.py (script to show the last image of the fractal deposit generated)
