set terminal post eps enhanced color "Helvetica" 18
set output "parplate1D_planewave_xyz_td.eps"
set title "Vulture Test Case: Parallel plate with partial planewave source: k_x, E_y, H_z"
set xlabel "Time (ns)"
set ylabel "Electric field, E_y(3,0,0) (V/m)"
plot "eh_op1_td.asc"                                                     us ($2/1e-9):4 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_planewave_xyz/eh_op1_td.asc" us ($2/1e-9):4 ti "Validation" w l ls 2
