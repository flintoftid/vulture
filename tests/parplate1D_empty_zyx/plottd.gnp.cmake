set terminal post eps enhanced color "Helvetica" 18
set output "parplate1D_empty_zyx_td.eps"
set title "Vulture Test Case: Parallel plate waveguide: k_z, E_y, H_x"
set xlabel "Time (ns)"
set ylabel "Electric field, E_y(0,0,100) (V/m)"
plot "eh_op1_td.asc"                                                 us ($2/1e-9):4 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_empty_zyx/eh_op1_td.asc" us ($2/1e-9):4 ti "Validation" w l ls 2
