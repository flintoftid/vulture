set terminal post eps enhanced color "Helvetica" 18
set output "softsrc_idy_fd.eps"
set title "Vulture Test Case: Soft source: IDY"
set xlabel "Frequency (MHz)"
set ylabel "Electric field, |E_y|(25,15.5,15) (dB V/m)"
plot "eh_op2_fd.asc"                                        us ($1/1e6):(10*log10($4**2+$5**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/softsrc_idy/eh_op2_fd.asc" us ($1/1e6):(10*log10($4**2+$5**2)) ti "Validation" w l ls 2

