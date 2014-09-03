set terminal post eps enhanced color "Helvetica" 18
set output "planewave_empty_xyz_td.eps"
set title "Vulture Test Case: Plane-wave in free-space: k_x, E_y, H_z"
set xlabel "Time (ns)"
set ylabel "Electric field, E_y(10,11,12) (V/m)"
plot "eh_op1_td.asc"                                                us ($2/1e-9):4 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/planewave_empty_xyz/eh_op1_td.asc" us ($2/1e-9):4 ti "Validation" w l ls 2
