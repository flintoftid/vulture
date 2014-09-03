set terminal post eps enhanced color "Helvetica" 18
set output "parplate1D_diel_td.eps"
set title "Vulture Test Case: Dielectric slab in parallel-plate waveguide"
set xlabel "Time (ns)"
set ylabel "Electric field, E_y(0,0,16) (V/m)"
plot "eh_op1_td.asc"                                            us ($2/1e-9):4 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_diel/eh_op1_td.asc" us ($2/1e-9):4 ti "Validation" w l ls 2
