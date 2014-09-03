set terminal post eps enhanced color "Helvetica" 18
set output "freespace_pml_fd.eps"
set title "Vulture Test Case: Freespace PML"
set xlabel "Frequency (MHz)"
set ylabel "Electric field, |E_z|(10,20,10) (dB V/m)"
plot "eh_op2_fd.asc"                                          us ($1/1e6):(10*log10($6**2+$7**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/freespace_pml/eh_op2_fd.asc" us ($1/1e6):(10*log10($6**2+$7**2)) ti "Validation" w l ls 2

