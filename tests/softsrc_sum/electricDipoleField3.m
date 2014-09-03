function [ E , H ] = electricDipoleField3( rvec , w , pvec )
%
% Field from electric dipole moment.
% Single space point, multiple frequencies.
%

  [ c0 , eps0 , mu0 , eta0 ] = emConst();

  assert( size( rvec ) == [ 1 ,3 ] );
  assert( size( w , 2 ) == 1 );
  assert( size( pvec , 2 ) == 3 );
  assert( size( w , 1 ) == size( pvec , 1 ) );

  % distance to observation point.
  r = sqrt( dot( rvec , rvec , 2 ) );

  % Unit vector along rvec.
  runit = repmat( rvec ./ r , [ size( pvec , 1 ) , 1 ] );

  % Frequency factors.
  k = w ./ c0;
  kr = k .* r;
  jkr = j .* kr;
  ejkrokr = exp( -jkr ) ./ kr;

  % Polarisation of H.
  rcrossp = cross( runit , pvec , 2 );

  % Polarisations of E.
  dir1 = cross( rcrossp , runit );
  dir2 = bsxfun( @times , 3 .* dot( runit , pvec , 2 ) , runit ) - pvec;

  % Fields.
  H = bsxfun( @times , c0 .* k.^3 ./ ( 4 .* pi ) .* ( 1 + 1 ./ jkr ) .* ejkrokr , rcrossp );
  E = bsxfun( @times , k.^3 ./ ( 4 .* pi .* eps0 ) .* ejkrokr , dir1 ) + ...
      bsxfun( @times , k.^3 ./ ( 4 .* pi .* eps0 ) .* ( 1 ./ kr.^2 + j ./ kr ) .* ejkrokr , dir2 );

end % function

