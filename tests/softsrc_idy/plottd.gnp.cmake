set terminal post eps enhanced color "Helvetica" 18
set output "softsrc_idy_td.eps"
set title "Vulture Test Case: Soft source: IDY"
set xlabel "Time (ns)"
set ylabel "Electric field, E_y(25,15.5,15) (V/m)"
plot "eh_op1_td.asc"                                        us ($2/1e-9):4 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/softsrc_idy/eh_op1_td.asc" us ($2/1e-9):4 ti "Validation" w l ls 2

