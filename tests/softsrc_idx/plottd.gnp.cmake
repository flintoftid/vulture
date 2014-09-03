set terminal post eps enhanced color "Helvetica" 18
set output "softsrc_idx_td.eps"
set title "Vulture Test Case: Soft source: IDX"
set xlabel "Time (ns)"
set ylabel "Electric field, E_x(15.5,15,25) (V/m)"
plot "eh_op1_td.asc"                                        us ($2/1e-9):3 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/softsrc_idx/eh_op1_td.asc" us ($2/1e-9):3 ti "Validation" w l ls 2

