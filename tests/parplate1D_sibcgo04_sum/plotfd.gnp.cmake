se terminal post eps enhanced color "Helvetica" 16
se output "parplate1D_sibcgo04_sum.eps"
se title "Vulture Test Case: Parallel plate WG GO04 SIBC: Summary"
se xlabel "Frequency (MHz)"
set ylabel "Magnitude of scattering parameter, |S_{ij}| (dB)"
set key top left
set log x
set xrange [1:10e3]
j=sqrt(-1)
plot "@VULTURE_SOURCE_DIR@/tests/parplate1D_sibcgo04_sum/GO04_Sfit.nport" us ($1/1e6):(db20(abs($10+j*$11)))  ti "Vector fit, S_{21}"    w l ls 1, \
     "@VULTURE_SOURCE_DIR@/tests/parplate1D_sibcgo04_sum/GO04_Sfit.nport" us ($1/1e6):(db20(abs($24+j*$25)))  ti "Vector fit, S_{34}"    w l ls 2, \
     "../parplate1D_sibcgo04_xyz/eh_op4_fd.asc"                           us ($1/1e6):(10*log10($4*$4+$5*$5)) ti "k_x, E_y, H_z: S_{21}" w l ls 3, \
     "../parplate1D_sibcgo04_xzy/eh_op4_fd.asc"                           us ($1/1e6):(10*log10($6*$6+$7*$7)) ti "k_x, E_z, H_y: S_{34}" w l ls 4, \
     "../parplate1D_sibcgo04_yxz/eh_op4_fd.asc"                           us ($1/1e6):(10*log10($2*$2+$3*$3)) ti "k_y, E_x, H_z: S_{34}" w l ls 5, \
     "../parplate1D_sibcgo04_yzx/eh_op4_fd.asc"                           us ($1/1e6):(10*log10($6*$6+$7*$7)) ti "k_y, E_z, H_x: S_{21}" w l ls 6, \
     "../parplate1D_sibcgo04_zxy/eh_op4_fd.asc"                           us ($1/1e6):(10*log10($2*$2+$3*$3)) ti "k_z, E_x, H_y: S_{21}" w l ls 7, \
     "../parplate1D_sibcgo04_zyx/eh_op4_fd.asc"                           us ($1/1e6):(10*log10($4*$4+$5*$5)) ti "k_z, E_y, H_x: S_{34}" w l ls 8

