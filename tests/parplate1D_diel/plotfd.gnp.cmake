set terminal post eps enhanced color "Helvetica" 18
set output "parplate1D_diel_fd.eps"
set title "Vulture Test Case: Dielectric slab in parallel-plate waveguide"
set xlabel "Frequency (MHz)"
set ylabel "Electric field, |E_y|(0,0,16) (dB V/m)"
set key bottom left
plot "eh_op2_fd.asc"                                            us ($1/1e6):(10*log10($4**2+$5**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_diel/eh_op2_fd.asc" us ($1/1e6):(10*log10($4**2+$5**2)) ti "Validation" w l ls 2, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_diel/analytic.dat"  us ($1/1e6):(20*log10($4))          ti "Analytic"   w l ls 3


