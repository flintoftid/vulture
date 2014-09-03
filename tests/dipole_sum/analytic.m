
f = linspace( 10e6 , 10e9 , 2000 )';

% Mesh size.
dl = 0.01;
len = ( 19 - 12 + 1 ) * dl;
diameter = 0.2 * dl;
Zl = 50.0;

[ AF , Zin , led , Rloss ] = afDipole( f , len , diameter , Zl );

writeDataFile( 'analytic.dat' , [ f , real( Zin ) , imag( Zin ) ] , { 'f [Hz]' , 'Re(Zin) [ohms]' , 'Im(Zin) (ohms)' } , { '' } );

