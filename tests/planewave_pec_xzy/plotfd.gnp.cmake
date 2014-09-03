set terminal post eps enhanced color "Helvetica" 18
set output "planewave_pec_xzy_fd.eps"
set title "Vulture Test Case: Plane-wave cut by PEC: k_x, E_z, H_y"
set xlabel "Frequency (MHz)"
set ylabel "Magnetic field, |H_y|(9,12,11) (dB A/m)"
plot "eh_op2_fd.asc"                                              us ($1/1e6):(10*log10($10**2+$11**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/planewave_pec_xzy/eh_op2_fd.asc" us ($1/1e6):(10*log10($10**2+$11**2)) ti "Validation" w l ls 2

