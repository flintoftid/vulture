db10( x ) = 10.0 * log10 ( abs( x )  )
db20( x ) = 20.0 * log10 ( abs( x )  )
db20ri( r , i ) = 10.0 * log10 ( r**2 + i**2 )
set terminal post eps enhanced color "Helvetica" 16
set output 'sgref_pml6_fd.eps'
set xlabel 'Frequency (GHz)'
set ylabel 'Local reflection error (dB)'
set key top left
set yrange [-160:0]
ref=-60
plot 'eh_xlo2_fd.asc' us ($1/1e9):(db20ri($4,$5)-ref) ti 'XLO, E_y' w lines, \
     'eh_xlo2_fd.asc' us ($1/1e9):(db20ri($6,$7)-ref) ti 'XLO, E_z' w lines, \
     'eh_xhi2_fd.asc' us ($1/1e9):(db20ri($4,$5)-ref) ti 'XHI, E_y' w lines, \
     'eh_xhi2_fd.asc' us ($1/1e9):(db20ri($6,$7)-ref) ti 'XHI, E_z' w lines, \
     'eh_ylo2_fd.asc' us ($1/1e9):(db20ri($2,$3)-ref) ti 'YLO, E_x' w lines, \
     'eh_ylo2_fd.asc' us ($1/1e9):(db20ri($6,$7)-ref) ti 'YLO, E_z' w lines, \
     'eh_yhi2_fd.asc' us ($1/1e9):(db20ri($2,$3)-ref) ti 'YHI, E_x' w lines, \
     'eh_yhi2_fd.asc' us ($1/1e9):(db20ri($6,$7)-ref) ti 'YHI, E_z' w lines, \
     'eh_zlo2_fd.asc' us ($1/1e9):(db20ri($2,$3)-ref) ti 'ZLO, E_x' w lines, \
     'eh_zlo2_fd.asc' us ($1/1e9):(db20ri($4,$5)-ref) ti 'ZLO, E_y' w lines, \
     'eh_zhi2_fd.asc' us ($1/1e9):(db20ri($2,$3)-ref) ti 'ZHI, E_x' w lines, \
     'eh_zhi2_fd.asc' us ($1/1e9):(db20ri($4,$5)-ref) ti 'ZHI, E_y' w lines

