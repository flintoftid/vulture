set terminal post eps enhanced color "Helvetica" 18
set output "planewave_empty_zxy_td.eps"
set title "Vulture Test Case: Plane-wave in free-space: k_z, E_x, H_y"
set xlabel "Time (ns)"
set ylabel "Electric field, E_x(11,12,10) (V/m)"
plot "eh_op1_td.asc"                                                us ($2/1e-9):3 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/planewave_empty_zxy/eh_op1_td.asc" us ($2/1e-9):3 ti "Validation" w l ls 2
