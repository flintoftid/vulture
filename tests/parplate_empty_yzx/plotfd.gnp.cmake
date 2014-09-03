set terminal post eps enhanced color "Helvetica" 18
set output "parplate_empty_yzx_fd.eps"
set title "Vulture Test Case: Parallel plate waveguide: k_y, E_z, H_x"
set xlabel "Frequency (MHz)"
set ylabel "Electric field, |E_z|(3,12,3) (dB V/m)"
plot "eh_op2_fd.asc"                                               us ($1/1e6):(10*log10($6**2+$7**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/parplate_empty_yzx/eh_op2_fd.asc" us ($1/1e6):(10*log10($6**2+$7**2)) ti "Validation" w l ls 2

