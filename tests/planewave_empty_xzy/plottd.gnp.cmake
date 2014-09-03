set terminal post eps enhanced color "Helvetica" 18
set output "planewave_empty_xzy_td.eps"
set title "Vulture Test Case: Plane-wave in free-space: k_x, E_z, H_y"
set xlabel "Time (ns)"
set ylabel "Electric field, E_z(10,12,11) (V/m)"
plot "eh_op1_td.asc"                                                us ($2/1e-9):5 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/planewave_empty_xzy/eh_op1_td.asc" us ($2/1e-9):5 ti "Validation" w l ls 2
