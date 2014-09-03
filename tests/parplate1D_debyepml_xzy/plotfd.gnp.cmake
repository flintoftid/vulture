set terminal post eps enhanced color "Helvetica" 18
set output "parplate1D_debyepml_xzy_fd.eps"
set title "Vulture Test Case: Parallel plate waveguide with Debye medium: k_x, E_z, H_y"
set xlabel "Frequency (MHz)"
set ylabel "Electric field, |E_z|(100,0,0) (dB V/m)"
plot "eh_op2_fd.asc"                                                    us ($1/1e6):(10*log10($6**2+$7**2)) ti "Test"       w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_debyepml_xzy/eh_op2_fd.asc" us ($1/1e6):(10*log10($6**2+$7**2)) ti "Validation" w l ls 2

