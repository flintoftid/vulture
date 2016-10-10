function [ A , f ] = tdfdFFT( a , wf , t )
%
%  [ A , f ] = tdfdFFT( a , t )
%
% FFT 3D field array.
% 
% Inputs:
%
% a           - field vector in time domain: a(i,j,k,ts) [-].
% wf          - reference waveform in time domain [-].
% t           - vector of sample times [s].
%
% Outputs:
%
% A           - field vector in frequency domain: A(i,j,k,f) [-].
% f           - vector of frequencies [Hz].
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

  Nt = length( t );
  dt = t(2) - t(1);

  df = 1 / ( Nt * dt );
  f = df.*(0:Nt-1).';

  A = fft( a , [] , 4 );

  WF(1,1,1,:) = fft( wf );

  Nf = Nt / 10;
  A = A(:,:,:,1:Nf) ./ repmat( WF(1,1,1,1:Nf) , [ size( A , 1) , size( A , 2 ) , size( A , 3 ) ] );
  f = f(1:Nf);
 
end % function

