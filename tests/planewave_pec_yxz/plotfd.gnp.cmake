set terminal post eps enhanced color "Helvetica" 18
set output "planewave_pec_yxz_fd.eps"
set title "Vulture Test Case: Plane-wave cut by PEC: k_y, E_x, H_z"
set xlabel "Frequency (MHz)"
set ylabel "Magnetic field, |H_z|(11,9,12) (dB A/m)"
plot "eh_op2_fd.asc"                                              us ($1/1e6):(10*log10($12**2+$13**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/planewave_pec_yxz/eh_op2_fd.asc" us ($1/1e6):(10*log10($12**2+$13**2)) ti "Validation" w l ls 2

