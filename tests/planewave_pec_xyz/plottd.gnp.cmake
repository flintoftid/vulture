set terminal post eps enhanced color "Helvetica" 18
set output "planewave_pec_xyz_td.eps"
set title "Vulture Test Case: Plane-wave cut by PEC: k_x, E_y, H_z"
set xlabel "Time (ns)"
set ylabel "MAgnetic field, H_z(9,11,12) (A/m)"
plot "eh_op1_td.asc"                                              us ($2/1e-9):8 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/planewave_pec_xyz/eh_op1_td.asc" us ($2/1e-9):8 ti "Validation" w l ls 2
