Written by Chaitanya Chintaluri, 01/12/2016
Updated on 08/12/2016

To replicate the results obtained for the FEM use Anaconda Scientific Python
distribution. For this, I provide my environment that can be re-created at your
desktop.

This file is provided as fenics27_spec.txt ; assumes that you use linux-64
conda create --name fenics27 --file fenics27_spec.txt
source activate fenics27

For the results shown in the plots in ./results/*.png the following was used.


1) Generate the mesh (gmsh 2.13.1 to create a .msh file from .geo, and dolfin
to create .xml files from .msh)

cd mesh
sh mesh_it_wm.sh
***NOW YOU WAIT***

2) Run FEM model

cd ..
python fem_4shell_wm.py
***NOW YOU WAIT***

(uses params from parameters_wm.py)
(gets you ./results/4Shell_FEM_*_wm.npy)

3) Analytical methods from Nunez, Appendix G

python analytical_4Shell_Nunez.py

(uses params from parameters_wm.py)
(gets you ./results/phi_*.npy)

4)Corrections by Solveig for Nunez Appendix G

python analytical_4Shell_Nunez_c.py

(uses params from parameters_wm.py)
(gets you ./results/phi_*_c.npy)

5) Make plots! 

source deactivate
python plot_potentials.py

(uses phi_*.npy, phi_*_c.npy, and 4Shell_FEM_*_wm.npy)
(gets you *.png)

NOTES:

1) fem_4shell.py; ./mesh/sphere_4.geo; ./mesh/mesh_it.geo
Can be used as before, but the meshes generated do not produce convergence to
the analytical solution for the 40 and 80 cases.

This is resolved by making another sphere inside the brain sphere, called wm
(analogous to white matter) but with the same conductivity as the brain. The
only difference is that the mesh sizes in wm are coarse - this resolves the
mis match with the analytical solution.

The resulting files are ./results/4Shell_FEM_*.npy (**no _wm**)

2) ./mesh/sphere_4_wm_high_res.geo (and change ./mesh_it_wm.geo and
fem_4shell_wm.py accordingly)

This mesh is used to establish that the numerical method has converged. By
increasing the fine-ness of the mess, and obtaining the results, the fine (high
resolution) are compared with the coarse resolution results 

These results are presented in ./results/4Shell_FEM_*_wm_high_res.npy

The differences between 4Shell_FEM_*_wm_high_res.npy and the corresponding
4Shell_FEM_*_wm.npy is ~ 1e-5.


