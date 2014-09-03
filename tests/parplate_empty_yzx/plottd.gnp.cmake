set terminal post eps enhanced color "Helvetica" 18
set output "parplate_empty_yzx_td.eps"
set title "Vulture Test Case: Parallel plate waveguide: k_y, E_z, H_x"
set xlabel "Time (ns)"
set ylabel "Electric field, E_z(3,12,3) (V/m)"
plot "eh_op4_td.asc"                                               us ($2/1e-9):5 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/parplate_empty_yzx/eh_op4_td.asc" us ($2/1e-9):5 ti "Validation" w l ls 2

