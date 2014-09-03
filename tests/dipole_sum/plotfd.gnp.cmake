set terminal post eps enhanced color "Helvetica" 18
set output "dipole_sum_fd.eps"
set title "Vulture Test Case: Dipole: Summary"
set xlabel "Frequency (GHz)"
set ylabel "Input impedance, Z_{in} (ohms)"
set yrange [-200:200]
set grid
dl=0.01
c0=299792458
dt=dl/2/c0
j = sqrt( -1 )
plot '@VULTURE_SOURCE_DIR@/tests/dipole_sum/analytic.dat' us ($1/1e9):($2)                                                   ti "Re(Z_{in}), analytic" w l ls 1, \
     '@VULTURE_SOURCE_DIR@/tests/dipole_sum/analytic.dat' us ($1/1e9):($3)                                                   ti "Im(Z_{in}), analytic" w l ls 2, \
     '../dipole_vrx/eh_op2_fd.asc'                        us ($1/1e9):(real(0.25*exp(j*2*pi*$1*dt/2)*($2+j*$3)/($10+j*$11))) ti "Re(Z_{in}), VRX"      w l ls 3, \
     '../dipole_vrx/eh_op2_fd.asc'                        us ($1/1e9):(imag(0.25*exp(j*2*pi*$1*dt/2)*($2+j*$3)/($10+j*$11))) ti "Im(Z_{in}), VRX"      w l ls 4, \
     '../dipole_vry/eh_op2_fd.asc'                        us ($1/1e9):(real(0.25*exp(j*2*pi*$1*dt/2)*($4+j*$5)/($12+j*$13))) ti "Re(Z_{in}), VRY"      w l ls 5, \
     '../dipole_vry/eh_op2_fd.asc'                        us ($1/1e9):(imag(0.25*exp(j*2*pi*$1*dt/2)*($4+j*$5)/($12+j*$13))) ti "Im(Z_{in}), VRY"      w l ls 6, \
     '../dipole_vrz/eh_op2_fd.asc'                        us ($1/1e9):(real(0.25*exp(j*2*pi*$1*dt/2)*($6+j*$7)/($8+j*$9)))   ti "Re(Z_{in}), VRZ"      w l ls 7, \
     '../dipole_vrz/eh_op2_fd.asc'                        us ($1/1e9):(imag(0.25*exp(j*2*pi*$1*dt/2)*($6+j*$7)/($8+j*$9)))   ti "Im(Z_{in}), VRZ"      w l ls 8

