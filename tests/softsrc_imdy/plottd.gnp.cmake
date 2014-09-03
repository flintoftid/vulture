set terminal post eps enhanced color "Helvetica" 18
set output "softsrc_imdy_td.eps"
set title "Vulture Test Case: Soft source: IMDY"
set xlabel "Time (ns)"
set ylabel "Electric field, E_x(25.5,15,15) (V/m)"
plot "eh_op1_td.asc"                                         us ($2/1e-9):3 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/softsrc_imdy/eh_op1_td.asc" us ($2/1e-9):3 ti "Validation" w l ls 2

