set terminal post eps enhanced color "Helvetica" 18
set output "parplate1D_sibcspar_xzy_fd.eps"
set title "Vulture Test Case: Parallel plate WG with S-parameter SIBC: k_x, E_z, H_y"
set xlabel "Frequency (MHz)"
set ylabel "Magnitude of scattering parameter, |S_{ij}| (dB)"
plot "eh_op2_fd.asc"                                                    us ($1/1e6):(10*log10($6**2+$7**2)) ti "Test, S_{11}"       w l ls 1, \
     "eh_op4_fd.asc"                                                    us ($1/1e6):(10*log10($6**2+$7**2)) ti "Test, S_{21}"       w l ls 2, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_sibcspar_xzy/eh_op2_fd.asc" us ($1/1e6):(10*log10($6**2+$7**2)) ti "Validation, S_{11}" w l ls 3, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_sibcspar_xzy/eh_op4_fd.asc" us ($1/1e6):(10*log10($6**2+$7**2)) ti "Validation, S_{21}" w l ls 4


