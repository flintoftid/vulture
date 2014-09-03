set terminal post eps enhanced color "Helvetica" 18
set output "planewave_pec_zxy_td.eps"
set title "Vulture Test Case: Plane-wave cut by PEC: k_z, E_x, H_y"
set xlabel "Time (ns)"
set ylabel "Magnetic field, H_y(11,12,9) (A/m)"
plot "eh_op1_td.asc"                                              us ($2/1e-9):7 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/planewave_pec_zxy/eh_op1_td.asc" us ($2/1e-9):7 ti "Validation" w l ls 2
