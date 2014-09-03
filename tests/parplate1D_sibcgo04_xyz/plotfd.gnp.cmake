set terminal post eps enhanced color "Helvetica" 18
set output "parplate1D_sibcgo04_xyz_fd.eps"
set title "Vulture Test Case: Parallel plate WG with GO04 SIBC: k_x, E_y, H_z"
set xlabel "Frequency (MHz)"
set ylabel "Magnitude of scattering parameter, |S_{ij}| (dB)"
plot "eh_op2_fd.asc"                                                    us ($1/1e6):(10*log10($4**2+$5**2)) ti "Test, S_{11}"       w l ls 1, \
     "eh_op4_fd.asc"                                                    us ($1/1e6):(10*log10($4**2+$5**2)) ti "Test, S_{21}"       w l ls 2, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_sibcgo04_xyz/eh_op2_fd.asc" us ($1/1e6):(10*log10($4**2+$5**2)) ti "Validation, S_{11}" w l ls 3, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_sibcgo04_xyz/eh_op4_fd.asc" us ($1/1e6):(10*log10($4**2+$5**2)) ti "Validation, S_{21}" w l ls 4

