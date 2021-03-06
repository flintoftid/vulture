set terminal post eps enhanced color "Helvetica" 18
set output "parplate1D_sibcgo04_xyz_td.eps"
set title "Vulture Test Case: Parallel plate WG with GO04 SIBC: k_x, E_y, H_z"
set xlabel "Time (ns)"
set ylabel "Electric field, E_y (V/m)"
plot "eh_op1_td.asc"                                                    us ($2/1e-9):4 ti "Test, R"       w l ls 1, \
     "eh_op3_td.asc"                                                    us ($2/1e-9):4 ti "Test, T"       w l ls 2, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_sibcgo04_xyz/eh_op1_td.asc" us ($2/1e-9):4 ti "Validation, R" w l ls 3, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_sibcgo04_xyz/eh_op3_td.asc" us ($2/1e-9):4 ti "Validation, T" w l ls 4
