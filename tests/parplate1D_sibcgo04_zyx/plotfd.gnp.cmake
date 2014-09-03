set terminal post eps enhanced color "Helvetica" 18
set output "parplate1D_sibcgo04_zyx_fd.eps"
set title "Vulture Test Case: Parallel plate WG with GO04 SIBC: k_z, E_y, H_x"
set xlabel "Frequency (MHz)"
set ylabel "Magnitude of scattering parameter, |S_{ij}| (dB)"
plot "eh_op2_fd.asc"                                                    us ($1/1e6):(10*log10($4**2+$5**2)) ti "Test, S_{11}"       w l ls 1, \
     "eh_op4_fd.asc"                                                    us ($1/1e6):(10*log10($4**2+$5**2)) ti "Test, S_{21}"       w l ls 2, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_sibcgo04_zyx/eh_op2_fd.asc" us ($1/1e6):(10*log10($4**2+$5**2)) ti "Validation, S_{11}" w l ls 3, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_sibcgo04_zyx/eh_op4_fd.asc" us ($1/1e6):(10*log10($4**2+$5**2)) ti "Validation, S_{21}" w l ls 4

