set terminal post eps enhanced color "Helvetica" 18
set output "planewave_pec_xyz_fd.eps"
set title "Vulture Test Case: Plane-wave cut by PEC: k_x, E_y, H_z"
set xlabel "Frequency (MHz)"
set ylabel "Magnetic field, |H_z|(9,11,12) (dB A/m)"
plot "eh_op2_fd.asc"                                              us ($1/1e6):(10*log10($12**2+$13**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/planewave_pec_xyz/eh_op2_fd.asc" us ($1/1e6):(10*log10($12**2+$13**2)) ti "Validation" w l ls 2

