set terminal post eps enhanced color "Helvetica" 18
set output "dipole_vrz_td.eps"
set title "Vulture Test Case: Dipole: VRZ"
set xlabel "Time (ns)"
set ylabel "Electric field, E_z(25,15,15.5) (V/m)"
plot "eh_op3_td.asc"                                       us ($2/1e-9):5 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/dipole_vrz/eh_op3_td.asc" us ($2/1e-9):5 ti "Validation" w l ls 2

