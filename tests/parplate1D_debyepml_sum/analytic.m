
[ c0 , eps0 , mu0 , eta0 ] = emConst();

debye=@(w,eps_r,sigma,r,p) eps_r + sigma ./ ( j .* w .* eps0 ) + r ./ ( j .* w - p )  + conj( r ) ./ ( j .* w - conj( p ) );

eps_r = 2.0;
mu_r = 1.0;
sigma = 0.01;
tau = 1e-9;
delEps = 2.0;
p = -1 / tau;
r = 0.5 * delEps / tau;

eta_L = eta0;
f = linspace( 1e6 , 6e9 , 1000 );
w = 2 .* pi .* f;
epsc_r = debye( w , eps_r , sigma , r , p );
eta_R = eta0 .* sqrt( mu_r ./ epsc_r );
tau_LR = 2 .* eta_R ./ ( eta_R + eta_L ); 
data = [ f' , abs(tau_LR)' ];
writeDataFile( 'analytic.dat' , data , { 'f [Hz]' , '|T| [-]' } , { 'Transmission coefficent of Debye medium' } );

