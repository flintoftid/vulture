set terminal post eps enhanced color "Helvetica" 18
set output "planewave_pec_yzx_fd.eps"
set title "Vulture Test Case: Plane-wave cut by PEC: k_y, E_z, H_x"
set xlabel "Frequency (MHz)"
set ylabel "Magnetic field, |H_x|(12,9,11) (dB A/m)"
plot "eh_op2_fd.asc"                                              us ($1/1e6):(10*log10($8**2+$9**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/planewave_pec_yzx/eh_op2_fd.asc" us ($1/1e6):(10*log10($8**2+$9**2)) ti "Validation" w l ls 2

