set terminal post eps enhanced color "Helvetica" 18
set output "parplate_empty_xzy_td.eps"
set title "Vulture Test Case: Parallel plate waveguide: k_x, E_z, H_y"
set xlabel "Time (ns)"
set ylabel "Electric field, E_z(3,12,3) (V/m)"
plot "eh_op4_td.asc"                                               us ($2/1e-9):5 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/parplate_empty_xzy/eh_op4_td.asc" us ($2/1e-9):5 ti "Validation" w l ls 2
