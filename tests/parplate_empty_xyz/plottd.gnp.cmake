set terminal post eps enhanced color "Helvetica" 18
set output "parplate_empty_xyz_td.eps"
set title "Vulture Test Case: Parallel plate waveguide: k_x, E_y, H_z"
set xlabel "Time (ns)"
set ylabel "Electric field, E_y(3,12,3) (V/m)"
plot "eh_op4_td.asc"                                               us ($2/1e-9):4 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/parplate_empty_xyz/eh_op4_td.asc" us ($2/1e-9):4 ti "Validation" w l ls 2

