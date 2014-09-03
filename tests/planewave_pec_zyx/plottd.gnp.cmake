set terminal post eps enhanced color "Helvetica" 18
set output "planewave_pec_zyx_td.eps"
set title "Vulture Test Case: Plane-wave cut by PEC: k_z, E_y, H_x"
set xlabel "Time (ns)"
set ylabel "Magnetic field, H_x(12,11,9) (A/m)"
plot "eh_op1_td.asc"                                              us ($2/1e-9):6 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/planewave_pec_zyx/eh_op1_td.asc" us ($2/1e-9):6 ti "Validation" w l ls 2
