set terminal post eps enhanced color "Helvetica" 18
set output "parplate1D_sibcgo04_yzx_fd.eps"
set title "Vulture Test Case: Parallel plate WG with GO04 SIBC: k_y, E_z, H_x"
set xlabel "Frequency (MHz)"
set ylabel "Magnitude of scattering parameter, |S_{ij}| (dB)"
plot "eh_op2_fd.asc"                                                    us ($1/1e6):(10*log10($6**2+$7**2)) ti "Test, S_{11}"       w l ls 1, \
     "eh_op4_fd.asc"                                                    us ($1/1e6):(10*log10($6**2+$7**2)) ti "Test, S_{21}"       w l ls 2, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_sibcgo04_yzx/eh_op2_fd.asc" us ($1/1e6):(10*log10($6**2+$7**2)) ti "Validation, S_{11}" w l ls 3, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_sibcgo04_yzx/eh_op4_fd.asc" us ($1/1e6):(10*log10($6**2+$7**2)) ti "Validation, S_{21}" w l ls 4

