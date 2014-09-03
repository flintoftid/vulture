set terminal post eps enhanced color "Helvetica" 18
set output "closedbox_go04_td.eps"
set title "Vulture Test Case: Closed GO04 box"
set xlabel "Time (ns)"
set ylabel "Electric field, E_x(8,8,8) (V/m)"
plot "eh_op1_td.asc"                                           us ($2/1e-9):3 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/closedbox_go04/eh_op1_td.asc" us ($2/1e-9):3 ti "Validation" w l ls 2

