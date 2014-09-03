set terminal post eps enhanced color "Helvetica" 18
set output "softsrc_imdz_td.eps"
set title "Vulture Test Case: Soft source: IMDZ"
set xlabel "Time (ns)"
set ylabel "Electric field, E_y(15,25.5,15) (V/m)"
plot "eh_op1_td.asc"                                         us ($2/1e-9):4 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/softsrc_imdz/eh_op1_td.asc" us ($2/1e-9):4 ti "Validation" w l ls 2

