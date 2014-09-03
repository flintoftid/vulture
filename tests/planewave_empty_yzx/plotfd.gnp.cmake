set terminal post eps enhanced color "Helvetica" 18
set output "planewave_empty_yzx_fd.eps"
set title "Vulture Test Case: Plane-wave in free-space: k_y, E_z, H_x"
set xlabel "Frequency (MHz)"
set ylabel "Electric field, |E_z|(12,10,11) (dB V/m)"
plot "eh_op2_fd.asc"                                                us ($1/1e6):(10*log10($6**2+$7**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/planewave_empty_yzx/eh_op2_fd.asc" us ($1/1e6):(10*log10($6**2+$7**2)) ti "Validation" w l ls 2

