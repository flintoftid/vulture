set terminal post eps enhanced color "Helvetica" 18
set output "planewave_pec_xzy_td.eps"
set title "Vulture Test Case: Plane-wave cut by PEC: k_x, E_z, H_y"
set xlabel "Time (ns)"
set ylabel "Magnetic field, H_y(9,12,11) (A/m)"
plot "eh_op1_td.asc"                                              us ($2/1e-9):7 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/planewave_pec_xzy/eh_op1_td.asc" us ($2/1e-9):7 ti "Validation" w l ls 2
