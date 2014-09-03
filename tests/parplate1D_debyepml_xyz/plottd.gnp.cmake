set terminal post eps enhanced color "Helvetica" 18
set output "parplate1D_debyepml_xyz_td.eps"
set title "Vulture Test Case: Parallel plate waveguide with Debye medium: k_x, E_y, H_z"
set xlabel "Time (ns)"
set ylabel "Electric field, E_y(100,0,0) (V/m)"
plot "eh_op1_td.asc"                                                    us ($2/1e-9):4 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_debyepml_xyz/eh_op1_td.asc" us ($2/1e-9):4 ti "Validation" w l ls 2
