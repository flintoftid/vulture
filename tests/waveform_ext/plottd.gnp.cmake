set terminal post eps enhanced color "Helvetica" 18
set output "freespace_mur_td.eps"
set title "Vulture Test Case: External waveform"
set xlabel "Time (ns)"
set ylabel "Electric field, E_z(10,20,10) (V/m)"
plot "eh_op1_td.asc"                                         us ($2/1e-9):5 ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/waveform_ext/eh_op1_td.asc" us ($2/1e-9):5 ti "Validation" w l ls 2

