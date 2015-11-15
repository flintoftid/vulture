#!/bin/sh

# Create picture of mesh.
gvulture -p enclosure_se.mesh
gnuplot mesh.gnp
mv mesh.eps enclosure_se-mesh.eps

# Create Gmsh mesh.
gvulture -p -m enclosure_se.mesh
mv mesh.msh enclosure_se.msh

# Run solver.
vulture -v enclosure_se.mesh

# Create picture of mesh with LS22 cube.
gvulture -p enclosure_se_ls22.mesh
gnuplot mesh.gnp
mv mesh.eps enclosure_se_ls22-mesh.eps

# Create Gmsh mesh.
gvulture -p -m enclosure_se_ls22.mesh
mv mesh.msh enclosure_se_ls22.msh

# Run solver.
vulture -v enclosure_se_ls22.mesh

# Plot results.
gnuplot plot_td.gnp
gnuplot plot_tdlog.gnp
gnuplot plot_fd.gnp

