set terminal post eps enhanced color "Helvetica" 18
set output "dipole_vrx_td.eps"
set title "Vulture Test Case: Dipole: VRX"
set xlabel "Time (ns)"
set ylabel "Electric field, E_x(15.5,25,15) (V/m)"
plot "eh_op3_td.asc"                                       us ($2/1e-9):5 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/dipole_vrx/eh_op3_td.asc" us ($2/1e-9):5 ti "Validation" w l ls 2

