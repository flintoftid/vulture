set terminal post eps enhanced color "Helvetica" 18
set output "softsrc_imdy_fd.eps"
set title "Vulture Test Case: Soft source: IMDY"
set xlabel "Frequency (MHz)"
set ylabel "Electric field, |E_x|(25.5,15,15) (dB V/m)"
plot "eh_op2_fd.asc"                                         us ($1/1e6):(10*log10($2**2+$3**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/softsrc_imdy/eh_op2_fd.asc" us ($1/1e6):(10*log10($2**2+$3**2)) ti "Validation" w l ls 2

