# High resolution mesh by default 
# For low res change sphere_4.geo to sphere_4_lowres.geo and sphere_4.msh to sphere_4_lowres.msh
gmsh -3 -optimize_netgen sphere_4.geo
dolfin-convert sphere_4.msh sphere_4.xml
