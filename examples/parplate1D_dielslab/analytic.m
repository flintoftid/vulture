
[ c0 , eps0 , mu0 , eta0 ] = emConst();

dl = 0.01;

eps_r = 2.0;
mu_r = 1.0;
sigma = 0.01;

f = linspace( 1e6 , 6e9 , 1000 );
w = 2 * pi .* f;

[ S ] = emMultiRef( f , eta0 , eta0 , [eps_r] , [sigma] , [mu_r] , [0.0] , [101*dl] , 'S' );

S11 = squeeze( S(1,1,:) );
S21 = squeeze( S(2,1,:) );

data = [ f(:) , abs( S11 ) , angled( S11 ) , abs( S21 ) , angled( S21 ) ]; 
writeDataFile( 'analytic.dat' , data , { 'f [Hz]' , '|S11| [-]' , '/_S11 (deg.)' , '|S21| (-)' , '/_S21 [deg.]' } , { 'Analytic solution' } );

