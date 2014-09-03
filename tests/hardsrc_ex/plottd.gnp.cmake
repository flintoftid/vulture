set terminal post eps enhanced color "Helvetica" 18
set output "hardsrc_ex_td.eps"
set title "Vulture Test Case: Hard source: EX"
set xlabel "Time (ns)"
set ylabel "Electric field, E_x(15.5,15,15,25) (V/m)"
plot "eh_op1_td.asc"                                       us ($2/1e-9):3 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/hardsrc_ex/eh_op1_td.asc" us ($2/1e-9):3 ti "Validation" w l ls 2

