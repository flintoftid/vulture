set terminal post eps enhanced color "Helvetica" 18
set output "softsrc_idz_td.eps"
set title "Vulture Test Case: Soft source: IDZ"
set xlabel "Time (ns)"
set ylabel "Electric field, E_z(15,25,15.5) (V/m)"
plot "eh_op1_td.asc"                                        us ($2/1e-9):5 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/softsrc_idz/eh_op1_td.asc" us ($2/1e-9):5 ti "Validation" w l ls 2

