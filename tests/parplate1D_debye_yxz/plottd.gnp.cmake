set terminal post eps enhanced color "Helvetica" 18
set output "parplate1D_debye_yxz_td.eps"
set title "Vulture Test Case: Parallel plate waveguide with Debye medium: k_y, E_x, H_z"
set xlabel "Time (ns)"
set ylabel "Electric field, E_x(0,100,0) (V/m)"
plot "eh_op1_td.asc"                                                 us ($2/1e-9):3 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_debye_yxz/eh_op1_td.asc" us ($2/1e-9):3 ti "Validation" w l ls 2
