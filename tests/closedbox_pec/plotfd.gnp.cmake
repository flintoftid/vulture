set terminal post eps enhanced color "Helvetica" 18
set output "closedbox_pec_fd.eps"
set title "Vulture Test Case: Closed PEC box"
set xlabel "Frequency (MHz)"
set ylabel "Electric field, |E_x|(7,7,7) (dB V/m)"
plot "eh_op2_fd.asc"                                          us ($1/1e6):(10*log10($2**2+$3**2)) ti "Test"       w l ls 1 , \
     "@VULTURE_SOURCE_DIR@/tests/closedbox_pec/eh_op2_fd.asc" us ($1/1e6):(10*log10($2**2+$3**2)) ti "Validation" w l ls 2

