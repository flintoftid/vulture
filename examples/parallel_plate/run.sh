#!/bin/sh

# Plot mesh.
gvulture -p parallel_plate.mesh
gnuplot mesh.gnp
mv mesh.eps parallel_plate_mesh.eps

# Run solver.
vulture -v parallel_plate.mesh

# Plot time response at mesh centre using ASCII observer.
gnuplot plot_td.gnp

# Plot frequency spectrum at mesh centre using ASCII observer.
gnuplot plot_fd.gnp

# Extract time series atr mesh centre from binary observer.
cat > process.dat << EOF
CE  Vulture example: An empty parallel-plate waveguide.
1
 0 20 2 0 40 2 10 10 2
1 400
0 5.99585e+08 1.49896e+06
1
 10 10 20 20 10 10 3
EOF

xtime

# Compare time series to ASCII observer.
gnuplot plot_td2.gnp

# Calculate frequency spectra, including phase.
xtransall phase

# Extract frequency spectrum at mesh centre from binary observer.
xfreq phase

# Compare to spectrum ASCII observer.
gnuplot plot_fd2.gnp

# Extract surface plot of time response at time-step 60.
cat > process.dat << EOF
CE  Vulture example: An empty parallel-plate waveguide.
1
 0 20 2 0 40 2 10 10 2
60 60
0 5.99585e+08 1.49896e+06
1
 0 20 0 40 10 10 3
EOF

xplane

# Plot time response over plane.
gnuplot plot_tds.gnp

