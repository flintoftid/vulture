set terminal post eps enhanced color "Helvetica" 16
set output "softsrcp_sum.eps"
set title "Vulture Test Case: Soft source: Electric summary"
set xlabel "Frequency (GHz)"
set ylabel "Electric field (dB V/m)"
set log x
plot "@VULTURE_SOURCE_DIR@/tests/softsrc_sum/analyticp.dat"  us ($1/1e9):(20*log10(abs($2)))     ti "Analytic" w l ls 1, \
     "../softsrc_idx/eh_op2_fd.asc"                          us ($1/1e9):(10*log10($2*$2+$3*$3)) ti "IDX"      w l ls 2, \
     "../softsrc_idy/eh_op2_fd.asc"                          us ($1/1e9):(10*log10($4*$4+$5*$5)) ti "IDY"      w l ls 2, \
     "../softsrc_idz/eh_op2_fd.asc"                          us ($1/1e9):(10*log10($6*$6+$7*$7)) ti "IDZ"      w l ls 4

set output "softsrcm_sum.eps"
set title "Vulture Test Case: Soft source: Magnetic summary"
set xlabel "Frequency (GHz)"
set ylabel "Electric field (dB V/m)"
set log x
plot "@VULTURE_SOURCE_DIR@/tests/softsrc_sum/analyticm.dat"  us ($1/1e9):(20*log10(abs($2)))     ti "Analytic" w l ls 1, \
     "../softsrc_imdx/eh_op2_fd.asc"                         us ($1/1e9):(10*log10($6*$6+$7*$7)) ti "IMDX"     w l ls 2, \
     "../softsrc_imdy/eh_op2_fd.asc"                         us ($1/1e9):(10*log10($2*$2+$3*$3)) ti "IMDY"     w l ls 3, \
     "../softsrc_imdz/eh_op2_fd.asc"                         us ($1/1e9):(10*log10($4*$4+$5*$5)) ti "IMDZ"     w l ls 4

