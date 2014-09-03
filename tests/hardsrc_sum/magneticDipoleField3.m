function [ E , H ] = magneticDipoleField3( rvec , w , mvec )
%
% Field from magnetic dipole moment.
% Single space point, multiple frequencies.
%

  [ c0 , eps0 , mu0 , eta0 ] = emConst();

  assert( size( rvec ) == [ 1 ,3 ] );
  assert( size( w , 2 ) == 1 );
  assert( size( mvec , 2 ) == 3 );
  assert( size( w , 1 ) == size( mvec , 1 ) );

  % distance to observation point.
  r = sqrt( dot( rvec , rvec , 2 ) );

  % Unit vector along rvec.
  runit = repmat( rvec ./ r , [ size( mvec , 1 ) , 1 ] );

  % Frequency factors.
  k = w ./ c0;
  kr = k .* r;
  jkr = j .* kr;
  ejkrokr = exp( -jkr ) ./ kr;

  % Polarisation of E.
  rcrossm = cross( runit , mvec , 2 );

  % Polarisations of H.
  dir1 = cross( rcrossm , runit );
  dir2 = bsxfun( @times , 3 .* dot( runit , mvec , 2 ) , runit ) - mvec;

  % Fields.
  E = bsxfun( @times , -eta0 .* k.^3 ./ ( 4 .* pi ) .* ( 1 + 1 ./ jkr ) .* ejkrokr , rcrossm );
  H = bsxfun( @times , k.^3 ./ ( 4 .* pi ) .* ejkrokr , dir1 ) + ...
      bsxfun( @times , k.^3 ./ ( 4 .* pi ) .* ( 1 ./ kr.^2 + j ./ kr ) .* ejkrokr , dir2 );

end % function

