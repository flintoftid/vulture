
c0 = 299792458;
mu0 = 4 * pi * 1e-7;
eps0 = 1.0 / ( mu0 * c0 * c0 );
eta0 = sqrt( mu0 / eps0 );

dl = 12e-5;
eps_r = 3.0;
mu_r = 1.0;
sigma = 0.0;
f = linspace(1e3, 5e11,2000);

[ S ] = emMultiRef( f , eta0 , eta0 , [eps_r] , [sigma] , [mu_r] , [0.0] , [(13-4+1)*dl] , 'S' );

S11 = squeeze( S(1,1,:) );
S21 = squeeze( S(2,1,:) );
data = [ f(:) , abs( S11 ) , angled( S11 ) , abs( S21 ) , angled( S21 ) ]; 
writeDataFile( 'analytic.dat' , data , { 'f [Hz]' , '|S11| [-]' , '/_S11 (deg.)' , '|S21| (-)' , '/_S21 [deg.]' } , { 'Example 4: Analytic solution' } );

