set terminal post eps enhanced color "Helvetica" 18
set output "parplate1D_simplepml_yzx_td.eps"
set title "Vulture Test Case: Parallel plate waveguide with simple medium: k_y, E_z, H_x"
set xlabel "Time (ns)"
set ylabel "Electric field, E_z(0,100,0) (V/m)"
plot "eh_op1_td.asc"                                                     us ($2/1e-9):5 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_simplepml_yzx/eh_op1_td.asc" us ($2/1e-9):5 ti "Validation" w l ls 2
