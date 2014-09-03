set terminal post eps enhanced color "Helvetica" 16
set output 'sgref_pml12_td.eps'
set xlabel 'Time (ns)'
set yrange [-160:0]
set key top left 
set ylabel 'Local reflection error (dB)'
ref=-60
plot 'eh_xlo1_td.asc' us ($2/1e-9):(db20($4)-ref) ti 'XLO, E_y' w lines, \
     'eh_xlo1_td.asc' us ($2/1e-9):(db20($5)-ref) ti 'XLO, E_z' w lines, \
     'eh_xhi1_td.asc' us ($2/1e-9):(db20($4)-ref) ti 'XHI, E_y' w lines, \
     'eh_xhi1_td.asc' us ($2/1e-9):(db20($5)-ref) ti 'XHI, E_z' w lines, \
     'eh_ylo1_td.asc' us ($2/1e-9):(db20($3)-ref) ti 'YLO, E_x' w lines, \
     'eh_ylo1_td.asc' us ($2/1e-9):(db20($5)-ref) ti 'YLO, E_z' w lines, \
     'eh_yhi1_td.asc' us ($2/1e-9):(db20($3)-ref) ti 'YHI, E_x' w lines, \
     'eh_yhi1_td.asc' us ($2/1e-9):(db20($5)-ref) ti 'YHI, E_z' w lines, \
     'eh_zlo1_td.asc' us ($2/1e-9):(db20($3)-ref) ti 'ZLO, E_x' w lines, \
     'eh_zlo1_td.asc' us ($2/1e-9):(db20($4)-ref) ti 'ZLO, E_y' w lines, \
     'eh_zhi1_td.asc' us ($2/1e-9):(db20($3)-ref) ti 'ZHI, E_x' w lines, \
     'eh_zhi1_td.asc' us ($2/1e-9):(db20($4)-ref) ti 'ZHI, E_y' w lines

