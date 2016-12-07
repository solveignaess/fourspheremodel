Written by Chaitanya Chintaluri, 01/12/2016

1) Generate the mesh (gmsh 2.13.1)

cd mesh
sh mesh_it.sh

***NOW YOU WAIT***

2) Setup fenics env in Anaconda as mine, use "fenics27_spec.txt" file ; assumes that you use linux-64 

conda create --name fenics27 --file fenics27_spec.txt
source activate fenics27

3) Run FEM model

python fem_4shell.py

***NOW YOU WAIT***

(uses params from parameters.py)
(gets you 4Shell_FEM_*.npy)

4) Analytical methods from Nunez, Appendix G

python analytical_4Shell_Nunez.py

(uses params from parameters.py)
(gets you phi_*.npy)

5)Corrections by Solveig for Nunez Appendix G

python analytical_4Shell_Nunez_c.py

(uses params from parameters.py)
(gets you phi_*_c.npy)

6) Make plots! 

python plot_potentials.py

(uses phi_*.npy, phi_*_c.npy, and 4Shell_FEM_*.npy)
(gets you *.png)

