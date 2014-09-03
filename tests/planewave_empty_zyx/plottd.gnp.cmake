set terminal post eps enhanced color "Helvetica" 18
set output "planewave_empty_zyx_td.eps"
set title "Vulture Test Case: Plane-wave in free-space: k_z, E_y, H_x"
set xlabel "Time (ns)"
set ylabel "Electric field, E_y(12,11,10) (V/m)"
plot "eh_op1_td.asc"                                                us ($2/1e-9):4 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/planewave_empty_zyx/eh_op1_td.asc" us ($2/1e-9):4 ti "Validation" w l ls 2
