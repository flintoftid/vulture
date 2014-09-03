set terminal post eps enhanced color "Helvetica" 18
set output "dipole_vry_td.eps"
set title "Vulture Test Case: Dipole: VRY"
set xlabel "Time (ns)"
set ylabel "Electric field, E_x(15,15.5,25) (V/m)"
plot "eh_op3_td.asc"                                       us ($2/1e-9):3 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/dipole_vry/eh_op3_td.asc" us ($2/1e-9):3 ti "Validation" w l ls 2

