Description
-----------

These are the scripts to obtain forward solution for the four sphere head
models to obtain the EEG potentials from the point dipole sources located
inside the brain of 4 different compositions (brain, csf, skull, scalp).

We implemented, the correct Analytical formulation for this, and an equivalent
numerical Finite element model, and provide with the previously presented
analytical methods. We validate our method Analytical and the numerical with
each other, and compare it with the previous ones. We also test the limiting
case when all the four spheres have the same conductivity.

Our presented results can be used for placing a dipole in any orientation, and
the electrodes at place (for analytical the electrodes must be > radial distance
of the dipole). These can be modified in the parameters.py file. While the
analytical solution is recommended for its fast computational time, the
numerical solution serves as a good starting point for further improvements for
a more sophisticated head models and more complicated conductivity profiles.

Authors
-------
Written by Chaitanya Chintaluri(1), Solveig Næss(2) and Torbjørn Ness(3)

1. Nencki Institute of Experimental Biology, Warsaw, Poland.
2. University of Oslo, Oslo, Norway.
3. Norwegian University of Life Sciences, Ås, Norway.

License
-------
Copyright 2017 Chaitanya Chintaluri, Solveig Næss, and Torbjørn Ness

This software is provided under the GNU General Purpose License version 3.0,
You will find a copy of this license within this folder, or from online here: 
https://www.gnu.org/licenses/gpl-3.0.txt


Packages
--------

To replicate the results obtained for the FEM use Anaconda Scientific Python
distribution. For this, I provide my environment that can be re-created at your
desktop. This file is provided as fenics27_spec.txt; assumes that you use linux-64

To obtain this environment. In your terminal type,  

conda create --name fenics27 --file fenics27_spec.txt
source activate fenics27

Read further about this from the Anaconda's managing environments webpage.

Additionally, you will need 
matplotlib python package
gmsh 2.13.1


Files and work flow
-------------------

1) Generate the mesh (USE gmsh 2.13.1 to create a .msh file from .geo, and dolfin
to create .xml files from .msh)

cd mesh
sh mesh_it_wm.sh
***NOW YOU WAIT***

2) Finite Element Method

cd ..
python numerical_fem.py
***NOW YOU WAIT***
(uses params from parameters.py)
(gets you ./results/Numerical_*.npz ; * is rad, tan, mix)

Computes the potentials on the scalp hemisphere for the cases of 
skull sigma = scalp sigma / 20; "fem_20"
skull sigma = scalp sigma / 40; "fem_40"
skull sigma = scalp sigma / 80; "fem_80"

For the radial, tangential and the 45 deg. oriented dipole

3) Analytical methods from Nunez and Srinivasan '06, Appendix G

python analytical_NunSri06.py
(uses params from parameters.py)
(gets you ./results/Analytical_NunSri06_rad.npz)

Computes the potentials on the scalp hemisphere for the cases of 
skull sigma = scalp sigma / 20; "phi_20"
skull sigma = scalp sigma / 40; "phi_40"
skull sigma = scalp sigma / 80; "phi_80"
All sigmas equal, case; "phi_lim"

For the radial dipole ONLY.

4) Analytical methods from Srinivasan '98, Appendix

python analytical_Sri98.py
(uses params from parameters.py)
(gets you ./results/Analytical_Sri98_rad.npz)

Computes the potentials on the scalp hemisphere for the cases of 
skull sigma = scalp sigma / 20; "phi_20"
skull sigma = scalp sigma / 40; "phi_40"
skull sigma = scalp sigma / 80; "phi_80"
All sigmas equal, case; "phi_lim"

For the radial dipole ONLY.

5) Analytical methods, Corrected solution

python analytical_correct.py !!!!***!!!!
(uses params from parameters.py)
(gets you ./results/Analytical_*.npz; *  is rad, tan, mix)

Computes the potentials on the scalp hemisphere for the cases of 
skull sigma = scalp sigma / 20; "fem_20"
skull sigma = scalp sigma / 40; "fem_40"
skull sigma = scalp sigma / 80; "fem_80"

For the radial, tangential and the 45 deg. oriented dipole

!!!!***!!!!
Alternate implementation available for this.

python CalcPotential4Sphere.py
(Self contained)


5) Make plots! 
source deactivate # matplotlib not available within the conda environment.

(i)

python figure_solveig_unit_upscaled.py
# Includes the scaling used in the manuscript

python figure2.py !!!!***!!!!
Compares that the analytical solution and the FEM solution have converged

!!!!***!!!!
An alternate view of this plot available see,
python figure2_alt.py

(ii)
python figure3.py

Compares the previous methods to the proposed solution.

(iii)
python figure4.py

Compares the limiting case of all equal conductivity for the analytical
methods.

(iv)
python plot_test.py

Compares the two analytical solution implementation for the various
configurations of the dipole orientations.

NOTES
-----

The mesh used ./mesh/sphere_4.geo actually has 5 spheres instead of 4.
The inner most sphere is there for the sake of simplifying the mesh size and
the convergence time for the FEM simulation.

A lower resolution of this mesh is also available ./mesh/sphere_4_lowres.geo
which establishes the convergence of the solution.
