set terminal post eps enhanced color "Helvetica" 18
set output "parplate1D_simplepml_yzx_fd.eps"
set title "Vulture Test Case: Parallel plate waveguide with simple medium: k_y, E_z, H_x"
set xlabel "Frequency (MHz)"
set ylabel "Electric field, |E_z|(0,100,0) (dB V/m)"
plot "eh_op2_fd.asc"                                                     us ($1/1e6):(10*log10($6**2+$7**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_simplepml_yzx/eh_op2_fd.asc" us ($1/1e6):(10*log10($6**2+$7**2)) ti "Validation" w l ls 2

