set terminal post eps enhanced color "Helvetica" 18
set output "planewave_empty_yzx_td.eps"
set title "Vulture Test Case: Plane-wave in free-space: k_y, E_z, H_x"
set xlabel "Time (ns)"
set ylabel "Electric field, E_z(12,10,11) (V/m)"
plot "eh_op1_td.asc"                                                us ($2/1e-9):5 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/planewave_empty_yzx/eh_op1_td.asc" us ($2/1e-9):5 ti "Validation" w l ls 2
