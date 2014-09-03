set terminal post eps enhanced color "Helvetica" 18
set output "hardsrc_ez_fd.eps"
set title "Vulture Test Case: Hard source: EZ"
set xlabel "Frequency (MHz)"
set ylabel "Electric field, |E_z|(15,25,15.5) (dB V/m)"
plot "eh_op2_fd.asc"                                       us ($1/1e6):(10*log10($6**2+$7**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/hardsrc_ez/eh_op2_fd.asc" us ($1/1e6):(10*log10($6**2+$7**2)) ti "Validation" w l ls 2

