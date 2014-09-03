function [ v ] = tdfdCompactPulse( t , width , delay )
%
%
%

%
% This file is part of Vulture.
%
% Vulture finite-difference time-domain electromagnetic solver.
% Copyright (C) 2011-2013 Ian David Flintoft
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

  dt = t(2) - t(1);

  if( nargin < 2 )
    width = 20 * dt;
  end %if

  if( nargin < 3 )
    delay = 0.0;
  end %if

  v = zeros( size( t ) );
  
  for k=1:length( t )
    if( ( t(k) - delay ) < 2.0 * width )
      wt = pi / width * ( t(k) - delay );
      v(k) = 1.0 / 32.0 * ( 10.0 - 15.0 * cos( wt ) + 6.0 * cos( 2.0 * wt ) - cos( 3.0 * wt ) );
    end %if
  end % for

end % function

