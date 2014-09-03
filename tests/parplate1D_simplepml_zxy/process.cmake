
file( WRITE  "process.dat" "CE  Vulture Test: Parallel plate waveguide with simple medium: k_z, E_x\n" )
file( APPEND "process.dat" "1\n" )
file( APPEND "process.dat" " 0 0 1 0 0 1 0 200 1\n" )
file( APPEND "process.dat" "1 5000\n" )
file( APPEND "process.dat" "0 5.99585e+10 1.19917e+07\n" )
file( APPEND "process.dat" "0.01\n" )
file( APPEND "process.dat" " 0 0 0 0 200 200 1\n" )

execute_process( COMMAND @XTIME_EXECUTABLE@ )
execute_process( COMMAND @XTRANSALL_EXECUTABLE@ phase )
execute_process( COMMAND @XFREQ_EXECUTABLE@ phase )

