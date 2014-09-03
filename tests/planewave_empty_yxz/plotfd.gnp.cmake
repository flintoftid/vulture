set terminal post eps enhanced color "Helvetica" 18
set output "planewave_empty_yxz_fd.eps"
set title "Vulture Test Case: Plane-wave in free-space: k_y, E_x, H_z"
set xlabel "Frequency (MHz)"
set ylabel "Electric field, |E_y|(11,10,12) (dB V/m)"
plot "eh_op2_fd.asc"                                                us ($1/1e6):(10*log10($2**2+$3**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/planewave_empty_yxz/eh_op2_fd.asc" us ($1/1e6):(10*log10($2**2+$3**2)) ti "Validation" w l ls 2

