
[ c0 , eps0 , mu0 , eta0 ] = emConst();

f = linspace( 100e6 , 10e9 , 2000 )';
w = 2 .* pi .* f;
jw = j .* w;

% Mesh size.
dl = 1e-2;

% Electric current moment.
Idz = ones( size( w ) );

% Electric dipole moment.
p_z = Idz ./ jw;
punit = [ 0 , 0 , 1 ];
evec = bsxfun( @times , p_z , punit ); 

% Location of dipole.
rvecp  = [ 15 , 15 , 15.5 ] .* dl;

% Observation point for E.
rvec = [ 15 , 25 , 15.5 ] .* dl;

[ E , H ]   = electricDipoleField3( rvec - rvecp  , w , evec );
Ez = E(:,3);

writeDataFile( 'analyticp.dat' , [ f , abs( Ez ) , angled( Ez ) ] , { 'f [Hz]' , '|E_z| [V/m]' , '/_ E_z(deg.)' } , { '' } );

% Magnetic current moment.
IMdz = ones( size( w ) );

% Magnetic dipole moment.
m_z = IMdz ./ jw ./ mu0;
munit = [ 0 , 0 , 1 ];
mvec = bsxfun( @times , m_z , munit );

% Location of dipole.
rvecm  = [ 15.5 , 15.5 , 15 ] .* dl;

% Observation point for E.
rvec = [ 15 , 25.5 , 15 ] .* dl;

[ E , H ] = magneticDipoleField3( rvec - rvecm  , w , mvec );
Ey = E(:,2);

writeDataFile( 'analyticm.dat' , [ f , abs( Ey ) , angled( Ey ) ] , { 'f [Hz]' , '|E_y| [V/m]' , '/_ E_y(deg.)' } , { '' } );

