set terminal post eps enhanced color "Helvetica" 18
set output "closedbox_pec_td.eps"
set title "Vulture Test Case: Closed PEC box"
set xlabel "Time (ns)"
set ylabel "Electric field, E_x(7,7,7) (V/m)"
plot "eh_op1_td.asc"                                          us ($2/1e-9):3 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/closedbox_pec/eh_op1_td.asc" us ($2/1e-9):3 ti "Validation" w l ls 2

