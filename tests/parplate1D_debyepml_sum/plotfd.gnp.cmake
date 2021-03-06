
se terminal post eps enhanced color "Helvetica" 16
se output "parplate1D_debyepml_sum.eps"
#se yra [-5:0]
se title "Vulture Test Case: Parallel plate waveguide with Debye medium: Summary"
se xlabel "Frequency (MHz)"
se ylabel "Electric field (dB V/m)"
se key bottom right

  plot "@VULTURE_SOURCE_DIR@/tests/parplate1D_debyepml_sum/analytic.dat" us   ($1/1e6):(20*log10($2)) ti "Analytic" w l, \
       "../parplate1D_debyepml_yzx/eh_op2_fd.asc" us ($1/1e6):(10*log10($6*$6+$7*$7)) ti "k_y, E_z, H_x: E_z(0,100,0)" w l, \
       "../parplate1D_debyepml_yxz/eh_op2_fd.asc" us ($1/1e6):(10*log10($2*$2+$3*$3)) ti "k_y, E_x, H_z: E_x(0,100,0)" w l, \
       "../parplate1D_debyepml_zxy/eh_op2_fd.asc" us ($1/1e6):(10*log10($2*$2+$3*$3)) ti "k_z, E_x, H_y: E_x(0,0,100)" w l, \
       "../parplate1D_debyepml_zyx/eh_op2_fd.asc" us ($1/1e6):(10*log10($4*$4+$5*$5)) ti "k_z, E_y, H_x: E_y(0,0,100)" w l, \
       "../parplate1D_debyepml_xyz/eh_op2_fd.asc" us ($1/1e6):(10*log10($4*$4+$5*$5)) ti "k_x, E_y, H_z: E_y(100,0,0)" w l, \
       "../parplate1D_debyepml_xzy/eh_op2_fd.asc" us ($1/1e6):(10*log10($6*$6+$7*$7)) ti "k_x, E_z, H_y: E_z(100,0,0)" w l 

se output

