set terminal post eps enhanced color "Helvetica" 18
set output "softsrc_idx_fd.eps"
set title "Vulture Test Case: Soft source: IDX"
set xlabel "Frequency (MHz)"
set ylabel "Electric field, |E_x|(15.5,15,25) (dB V/m)"
plot "eh_op2_fd.asc"                                        us ($1/1e6):(10*log10($2**2+$3**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/softsrc_idx/eh_op2_fd.asc" us ($1/1e6):(10*log10($2**2+$3**2)) ti "Validation" w l ls 2

