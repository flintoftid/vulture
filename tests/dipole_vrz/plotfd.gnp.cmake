set terminal post eps enhanced color "Helvetica" 18
set output "dipole_vrz_fd.eps"
set title "Vulture Test Case: Dipole: VRZ"
set xlabel "Frequency (GHz)"
set ylabel "Input impedance, Z_{in} (ohms)"
set yrange [-200:200]
set grid
dl=0.01
c0=299792458
dt=dl/2/c0
j = sqrt( -1 )
  plot 'eh_op2_fd.asc'                                       us ($1/1e9):(real(0.25*exp(j*2*pi*$1*dt/2)*($6+j*$7)/($8+j*$9))) ti "Re(Z_{in}), test"       w l ls 1, \
       'eh_op2_fd.asc'                                       us ($1/1e9):(imag(0.25*exp(j*2*pi*$1*dt/2)*($6+j*$7)/($8+j*$9))) ti "Im(Z_{in}), test"       w l ls 2, \
       '@VULTURE_SOURCE_DIR@/tests/dipole_vrz/eh_op2_fd.asc' us ($1/1e9):(real(0.25*exp(j*2*pi*$1*dt/2)*($6+j*$7)/($8+j*$9))) ti "Re(Z_{in}), validation" w l ls 3, \
       '@VULTURE_SOURCE_DIR@/tests/dipole_vrz/eh_op2_fd.asc' us ($1/1e9):(imag(0.25*exp(j*2*pi*$1*dt/2)*($6+j*$7)/($8+j*$9))) ti "Re(Z_{in}), validation" w l ls 4

