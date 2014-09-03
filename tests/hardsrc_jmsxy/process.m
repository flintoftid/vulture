
c0 = 299792458; 

opNum = 1;
freqNum = 10;

% Magnitude of electric field in aperture.
E0 = 1.0;

% FDTD mesh intervals.
dx = 0.01;
dy = 0.01;
dz = 0.01;

% Aperture side and diagonal lengths.
a = ( 20 - 10 + 1 ) * dx;
b = ( 18 - 12     ) * dy;
D = sqrt( a^2 + b^2 );

% Load frequency domain planes in obervation plane.
[ ii , jj , kk , Ex , Ey , Ez , Hx , Hy , Hz , f ] = tdfdReadFrequencyDat3D( opNum , freqNum , 1 );

% Determine Cartesian coordinates of Ey edges centred on aperture.
x = ( ii - 15 ) .* dx;
y = ( jj - 15 + 0.5 ) .* dy;
z = ( kk - 0 ) .* dz;

% Set up grid data for plotting.
[ xx , yy ] = meshgrid( x ,y );
zz = z(1) .* ones( size( xx ) );

% FDTD fields in decibels.
ExdB = 20 .* log10( abs( Ex ) );
EydB = 20 .* log10( abs( Ey ) );
EzdB = 20 .* log10( abs( Ez ) );
EtdB = 20 .* log10( abs( Ex ).^2 + abs( Ey ).^2 + abs( Ez ).^2 );

% Wavelength, angular frequency, wave number.
lambda = c0 / f;
w = 2 * pi * f;
k = w / c0;
jk = j * k;

% Rayleigh range.
R = 2 * D^2 / lambda;

% Spherical coordinates of observation points. 
rr = sqrt( xx.^2 + yy.^2 + zz.^2 );
theta = acos( zz ./ rr );
phi = atan2( yy , xx );

% Spherical components of wave vector.
X = 0.5 .* k .* a .* sin( theta ) .* cos( phi );
Y = 0.5 .* k .* b .* sin( theta ) .* sin( phi );
v_x = X ./ pi;
v_y = Y ./ pi;
f_y = sinc( v_x ) .* sinc( v_y );

% Spherical components of radiated field.
E_r = 0.0 .* rr;
E_theta = jk .* a .* b .* E0 .* exp( -jk .* rr ) ./ ( 2 .* pi .* rr ) .* sin( phi ) .* f_y;
E_phi   = jk .* a .* b .* E0 .* exp( -jk .* rr ) ./ ( 2 .* pi .* rr ) .* cos( theta ) .* cos( phi ) .* f_y;

% Cartesian components of radiated field.
Exa = E_r .* sin( theta ) .* cos( phi ) + E_theta .* cos( theta ) .* cos( phi ) - E_phi .* sin( phi );
Eya = E_r .* sin( theta ) .* sin( phi ) + E_theta .* cos( theta ) .* sin( phi ) + E_phi .* cos( phi );
Eza = E_r .* cos( theta ) - E_theta .* sin( theta );
Eta = sqrt ( abs( E_r ).^2 + abs( E_theta ).^2 + abs( E_phi ).^2 );

% Analytic field in decibels.
ExadB = 20 .* log10( abs( Exa ) );
EyadB = 20 .* log10( abs( Eya ) );
EzadB = 20 .* log10( abs( Eza ) );
EtadB = 20 .* log10( abs( Eta ) );

% Maximum values of FDTD and analytic field and their difference.
EydBmax = max( max( EydB ) )
EyadBmax = max( max( EyadB ) )
offset = EydBmax - EyadBmax

%
% Plots.
%

figure( 1 );
hl1 = surf( xx , yy , EydB );
hold on;
hl2 = surf( xx , yy , EyadB - 12 );
hxl = xlabel( 'x (m)' );
hyl = ylabel( 'y (m)' );
hzl = zlabel( 'E_y (dB V/m)' );
print( 'hardsrc_jmsxy_surf' , '-depsc2' );
hold off;

%
%
%

clear all;

c0 = 299792458;
mu0 = 4 * pi * 1e-7;
eps0 = 1.0 / ( mu0 * c0 * c0 );
eta0 = sqrt( mu0 / eps0 );

opNum = 1;
freqNum = 1:50;

% Magnitude of electric field in aperture.
E0 = 1.0;

% FDTD mesh intervals.
dx = 0.01;
dy = 0.01;
dz = 0.01;

% Aperture side and diagonal lengths.
a = ( 20 - 10 + 1 ) * dx;
b = ( 18 - 12 + 1 ) * dy;
D = sqrt( a^2 + b^2 );

% Load frequency domain planes in obervation plane.
[ ii , jj , kk , Ex , Ey , Ez , Hx , Hy , Hz , f ] = tdfdReadFrequencyDat3D( opNum , freqNum , 1 );
Ex = squeeze( Ex(15,15,1,:) );
Ey = squeeze( Ey(15,15,1,:) );
Ez = squeeze( Ez(15,15,1,:) );
ii = 15;
jj = 15;
f = f(:);

% Determine Cartesian coordinates of Ey edges centred on aperture.
x = ( ii - 15 ) .* dx;
y = ( jj - 15 + 0.5 ) .* dy;
z = ( kk - 0 ) .* dz;

% FDTD fields in decibels.
ExdB = 20 .* log10( abs( Ex ) );
EydB = 20 .* log10( abs( Ey ) );
EzdB = 20 .* log10( abs( Ez ) );
EtdB = 20 .* log10( abs( Ex ).^2 + abs( Ey ).^2 + abs( Ez ).^2 );

% Wavelength, angular frequency, wave number.
lambda = c0 ./ f;
w = 2 .* pi .* f;
k = w ./ c0;
jk = j .* k;

% Rayleigh range.
R = 2 .* D.^2 ./ lambda;

% Spherical coordinates of observation point. 
r = sqrt( x^2 + y^2 + z^2 );
theta = acos( z / r );
phi = atan2( y , x );

% Spherical components of wave vector.
X = 0.5 .* k .* a .* sin( theta ) .* cos( phi );
Y = 0.5 .* k .* b .* sin( theta ) .* sin( phi );
v_x = X ./ pi;
v_y = Y ./ pi;
f_y = sinc( v_x ) .* sinc( v_y );

% Spherical components of radiated field.
E_r = 0.0 .* r;
E_theta = jk .* a .* b .* E0 .* exp( -jk .* r ) ./ ( 2 .* pi .* r ) .* sin( phi ) .* f_y;
E_phi   = jk .* a .* b .* E0 .* exp( -jk .* r ) ./ ( 2 .* pi .* r ) .* cos( theta ) .* cos( phi ) .* f_y;

% Cartesian components of radiated field.
Exa = E_r .* sin( theta ) .* cos( phi ) + E_theta .* cos( theta ) .* cos( phi ) - E_phi .* sin( phi );
Eya = E_r .* sin( theta ) .* sin( phi ) + E_theta .* cos( theta ) .* sin( phi ) + E_phi .* cos( phi );
Eza = E_r .* cos( theta ) - E_theta .* sin( theta );
Eta = sqrt ( abs( E_r ).^2 + abs( E_theta ).^2 + abs( E_phi ).^2 );

% Analytic field in decibels.
ExadB = 20 .* log10( abs( Exa ) );
EyadB = 20 .* log10( abs( Eya ) );
EzadB = 20 .* log10( abs( Eza ) );
EtadB = 20 .* log10( abs( Eta ) );

%
% Plots.
%

figure( 2 );
hl1 = semilogx( f / 1e6 , EydB , 'r-o' );
hold on;
hl2 = semilogx( f / 1e6 , EyadB - 12 , 'b-^' );
hxl = xlabel( 'Frequency (MHz)' );
hyl = ylabel( 'E_y(15,15,25) (dB V/m)' );
legend( [ hl1 , hl2 ] , 'FDTD' , 'Analytic' );
axis( [ 10 , 10000 , -80 , 0 ] );
print( 'hardsrc_jmsxy_axis' , '-depsc2' );
hold off;

figure( 3 );
hl1 = loglog( f / 1e6 , r ./ R , 'r-o' );
hold on;
hl2 = loglog( f / 1e6 , ones( size( r ./ R ) )  , 'r:' );
hxl = xlabel( 'Frequency (MHz)' );
hyl = ylabel( 'Distance / Rayleigh range (-)' );
axis( [ 10 , 10000 ] );
print( 'hardsrc_jmsxy_rayrange' , '-depsc2' );
hold off;
