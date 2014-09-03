function [ v ] = tdfdGaussianPulse( t , width , delay )
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
    width = 5 * sqrt( 2 ) * dt;
  end %if

  if( nargin < 3 )
    delay = 40 * dt;
  end %if

  v = exp( -0.5 .* ( ( ( t - delay ) ./ width ).^2 ) );

end % function

