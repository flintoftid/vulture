se terminal post eps enhanced color "Helvetica" 16
se output "parplate1D_sibcspar_sum.eps"
se title "Vulture Test Case: Parallel plate WG with S-parameter SIBC: Summary"
se xlabel "Frequency (MHz)"
set ylabel "Magnitude of scattering parameter, |S_{ij}| (dB)"
se key top left
plot -6                                                                             ti "Analytic, S_{11}"      w l, \
     "../parplate1D_sibcspar_yzx/eh_op2_fd.asc" us ($1/1e6):(10*log10($6*$6+$7*$7)) ti "k_y, E_z, H_x: S_{11}" w l, \
     "../parplate1D_sibcspar_yxz/eh_op2_fd.asc" us ($1/1e6):(10*log10($2*$2+$3*$3)) ti "k_y, E_x, H_z: S_{11}" w l, \
     "../parplate1D_sibcspar_zxy/eh_op2_fd.asc" us ($1/1e6):(10*log10($2*$2+$3*$3)) ti "k_z, E_x, H_y: S_{11}" w l, \
     "../parplate1D_sibcspar_zyx/eh_op2_fd.asc" us ($1/1e6):(10*log10($4*$4+$5*$5)) ti "k_z, E_y, H_x: S_{11}" w l, \
     "../parplate1D_sibcspar_xyz/eh_op2_fd.asc" us ($1/1e6):(10*log10($4*$4+$5*$5)) ti "k_x, E_y, H_z: S_{11}" w l, \
     "../parplate1D_sibcspar_xzy/eh_op2_fd.asc" us ($1/1e6):(10*log10($6*$6+$7*$7)) ti "k_x, E_z, H_y: S_{11}" w l, \
     -12                                                                            ti "Analytic, S_{21}"      w l, \
     "../parplate1D_sibcspar_yzx/eh_op4_fd.asc" us ($1/1e6):(10*log10($6*$6+$7*$7)) ti "k_y, E_z, H_x: S_{21}" w l, \
     "../parplate1D_sibcspar_yxz/eh_op4_fd.asc" us ($1/1e6):(10*log10($2*$2+$3*$3)) ti "k_y, E_x, H_z: S_{21}" w l, \
     "../parplate1D_sibcspar_zxy/eh_op4_fd.asc" us ($1/1e6):(10*log10($2*$2+$3*$3)) ti "k_z, E_x, H_y: S_{21}" w l, \
     "../parplate1D_sibcspar_zyx/eh_op4_fd.asc" us ($1/1e6):(10*log10($4*$4+$5*$5)) ti "k_z, E_y, H_x: S_{21}" w l, \
     "../parplate1D_sibcspar_xyz/eh_op4_fd.asc" us ($1/1e6):(10*log10($4*$4+$5*$5)) ti "k_x, E_y, H_z: S_{21}" w l, \
     "../parplate1D_sibcspar_xzy/eh_op4_fd.asc" us ($1/1e6):(10*log10($6*$6+$7*$7)) ti "k_x, E_z, H_y: S_{21}" w l

