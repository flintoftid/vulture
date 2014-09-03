set terminal post eps enhanced color "Helvetica" 16
set output "hardsrcp_sum.eps"
set title "Vulture Test Case: Hard source: Electric summary"
set xlabel "Frequency (GHz)"
set ylabel "Electric field (dB V/m)"
set log x
plot "@VULTURE_SOURCE_DIR@/tests/hardsrc_sum/analyticp.dat"  us ($1/1e9):(20*log10(abs($2)))     ti "Analytic" w l ls 1, \
     "../hardsrc_ex/eh_op2_fd.asc"                           us ($1/1e9):(10*log10($2*$2+$3*$3)) ti "=EX"      w l ls 2, \
     "../hardsrc_ey/eh_op2_fd.asc"                           us ($1/1e9):(10*log10($4*$4+$5*$5)) ti "=EY"      w l ls 2, \
     "../hardsrc_ez/eh_op2_fd.asc"                           us ($1/1e9):(10*log10($6*$6+$7*$7)) ti "=EZ"      w l ls 4

