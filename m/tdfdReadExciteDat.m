function [ wf ] = tdfdReadExciteDat( t )
%
%
%

%
% This file is part of Vulture.
%
% Vulture finite-difference time-domain electromagnetic solver.
% Copyright (C) 2011-2016 Ian David Flintoft
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software Foundation,
% Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
%

% Author: I. D. Flintoft 
% Date: [FIXME]

  % Time-step.
  dt = t(2) - t(1);

  % Name of exciation file.
  exciteFileName = 'excite.dat';

  % Open file.
  fd = fopen( exciteFileName , 'r' );
  if( fd == -1 )
    error( 'Error: Cannot open file %s.' , exciteFileName );
  end %if

  % Get excitation type.
  [ excitationType , count ] = fscanf( fd , '%d' , 1 );
  if( count ~= 1 )
    error( 'Error reading from %s.' , exciteFileName );
  end % if

  % Construct waveform for each type.
  switch( excitationType ) 
  case 0 % IMPULSE
    wf = zeros( 1 , size( t ) );
    wf(1) = 1.0;
    break;
  case 1 % GAUSSIAN
  case 3 % SHORT_GAUSSIAN
    [ data , count ] = fscanf( fd , '%e' , 2 );
    alpha = data(1);
    beta = data(2);
    wf = exp( -( alpha .* ( t - beta ) ).^2 );
    break;
  case 2 %SINUSOIDAL
    [ data , count ] = fscanf( fd , '%e' , 1 );
    omega = 2.0 * pi * data(1);
    wf = sin( omega .* t );
    break;
  case 4 % CONSTANT
    [ data , count ] = fscanf( fd , '%e' , 2 );
    alpha = data(1);
    beta = data(2);    
    gaussint = 0.0;
    for ts=1:length( t )
      gaussint = gaussint + exp( -( alpha * ( t(ts) - beta ) )^2 );
      wf(ts) = gaussint * dt * alpha / sqrt( pi );
    end %for
    break;
  case 5 %GAUSSIAN_MODULATED_SINUSOID
    [ data , count ] = fscanf( fd , '%e' , 3 );
    alpha = data(1);
    beta = data(2);
    omega = 2.0 * pi * data(3);
    wf = exp( -( alpha .* ( t - beta ) ).^2 ) .* sin( omega .* t );
    break;
  case 6 % DIFFERENTIATED_GAUSSIAN
    [ data , count ] = fscanf( fd , '%e' , 2 );
    alpha = data(1);
    beta = data(2); 
    wf = -exp( -( alpha .* ( t - beta ) ).^2 ) .* 2.0 .* alpha .* ( t - beta );
    break;
  case 7 %ESD_PULSE
    tns = t * 1e9;
    wf = ( 7 .* ( exp( -tns ./ 1.5 ) - exp( -tns ./ 0.6 ) ) + ...
         1.6 .* ( exp( -tns ./ 26 ) - exp( -tns ./ 3.5 ) ) ) ./ 2.61;
    break;
  case 8 % USER_DEFINED
    [ wf , count ] =fscanf( fd , '%e' , size( t ) );
    if( count ~= length( t ) )
      error( 'Error reading user excitation data from %s.' , exciteFileName );
    end % if
    break;
  otherwise
    error( 'Error: excitation type number %d invalid.' , excitationType );
  end % switch

  % Close file.
  fclose( fd );

end % function
