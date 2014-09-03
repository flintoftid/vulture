set terminal post eps enhanced color "Helvetica" 18
set output "planewave_empty_xyz_fd.eps"
set title "Vulture Test Case: Plane-wave in free-space: k_x, E_y, H_z"
set xlabel "Frequency (MHz)"
set ylabel "Electric field, |E_y|(10,11,12) (dB V/m)"
plot "eh_op2_fd.asc"                                                us ($1/1e6):(10*log10($4**2+$5**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/planewave_empty_xyz/eh_op2_fd.asc" us ($1/1e6):(10*log10($4**2+$5**2)) ti "Validation" w l ls 2

