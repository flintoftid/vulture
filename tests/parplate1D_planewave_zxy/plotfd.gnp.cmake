set terminal post eps enhanced color "Helvetica" 18
set output "parplate1D_planewave_zxy_fd.eps"
set title "Vulture Test Case: Parallel plate with partial planewave source: k_z, E_x, H_y"
set xlabel "Frequency (MHz)"
set ylabel "Electric field, |E_x|(0,0,3) (dB V/m)"
plot "eh_op2_fd.asc"                                                     us ($1/1e6):(10*log10($2**2+$3**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_planewave_zxy/eh_op2_fd.asc" us ($1/1e6):(10*log10($2**2+$3**2)) ti "Validation" w l ls 2

