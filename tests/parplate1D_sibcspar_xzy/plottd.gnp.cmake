set terminal post eps enhanced color "Helvetica" 18
set output "parplate1D_sibcspar_xzy_td.eps"
set title "Vulture Test Case: Parallel plate WG with S-parameter SIBC: k_x, E_z, H_y"
set xlabel "Time (ns)"
set ylabel "Electric field, E_z (V/m)"
plot "eh_op1_td.asc"                                                    us ($2/1e-9):5 ti "Test, R"       w l ls 1, \
     "eh_op3_td.asc"                                                    us ($2/1e-9):5 ti "Test, T"       w l ls 2, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_sibcspar_xzy/eh_op1_td.asc" us ($2/1e-9):5 ti "Validation, R" w l ls 3, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_sibcspar_xzy/eh_op3_td.asc" us ($2/1e-9):5 ti "Validation, T" w l ls 4
