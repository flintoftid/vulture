
[ c0 , eps0 , mu0 , eta0 ] = emConst();

f = linspace( 10e6 , 1.6e9 , 500 );
w = 2 * pi * f;
k0 = w ./ c0;
lambda = c0 ./ f;

a = 300e-3;
b = 120e-3;
d = 300e-3;

l = 200e-3;
w = 30e-3;
t = 1.5e-3;

p = 150e-3;

v0 = 1;
Z0 = eta0;

we = w - ( 5 * t ) / ( 4 * pi ) * ( 1 + log( 4 * pi * w / t) );

Z0s = 120 * pi^2 / log( 2 * ( 1 + ( 1 - ( we / b )^2 )^0.25 ) / ( 1 - ( 1 - ( we / b )^2 )^0.25 ) );

Zap = 0.5 .* l ./ a .* j .* Z0s .* tan( k0 .* l ./ 2 );

v1 = v0 .* Zap ./ ( Z0 + Zap );
Z1 = Z0 .* Zap ./ ( Z0 + Zap );
Zg = Z0 ./ sqrt( 1 - ( lambda ./ 2 ./ a ).^2 );
kg = k0 .* sqrt( 1 - ( lambda ./ 2 ./ a ).^2 );
v2 = v1 ./ ( cos( kg .* p ) + j .* Z1 ./ Zg .* sin( kg .* p ) );
Z2 = ( Z1 + j .* Zg .* tan( kg .* p ) ) ./ ( 1 + j .* Z1 ./ Zg .* tan( kg .* p ) );
Z3 = j .* Zg .* tan( kg .* ( d - p ) );
vp = v2 .* Z3 ./ ( Z2 + Z3 );
ip = v2 ./ ( Z2 + Z3 );
vpp = v0 / 2;
ipp = v0 / 2 / Z0;

SE_E = abs( vpp ./ vp );
SE_H = abs( ipp ./ ip );

data = [ f' , SE_E' , SE_H' ];
writeDataFile( 'analytic.dat' , data , { 'f [Hz]' , 'SE_E [-]' , 'SE_H [-]' } , { 'Vulture example: Enclosure SE' } );
 

