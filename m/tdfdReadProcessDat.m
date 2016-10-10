function [ header , ko , io , nt1 , nt2 , d1 , d2 , ds , meshSize , ...
           ilo , ihi , jlo , jhi , klo , khi , fieldOrient ] = ...
    tdfdReadProcessDat( )
%
%  [ ko , io , nt1 , nt2 , d1 , d2 , ds , meshSize , ...
%    ilo , ihi , jlo , jhi , klo , khi , fieldOrient ] = ...
%                                       tdfdReadProcessDat( )
% 
%
% Read in processing control data from time domain code 
% process.dat file.
%
% Inputs:
% 
% None.
%
% Outputs:
%
% header - file header string
% ko - number of OP lines (-)
% io - array of cells for each OP lines:
%
%       io(n=1...ko,1)  xlo for OP line n (cells)
%       io(n=1...ko,2)  xhi for OP line n (cells)
%       io(n=1...ko,3)  xstep for OP line n (cells)
%       io(n=1...ko,4)  ylo for OP line n (cells)
%       io(n=1...ko,5)  yhi for OP line n (cells)
%       io(n=1...ko,6)  ystep for OP line n (cells)
%       io(n=1...ko,7)  zlo for OP line n (cells)
%       io(n=1...ko,8)  zhi for OP line n (cells)
%       io(n=1...ko,9)  zstep for OP line n (cells)
%
% nt1 - start frequenct step number (-)
% nt2 - stop frequency step number (-)
% d1 - start frequency (Hz)
% d2 - stop frequency (Hz)
% ds - frequency step (Hz)
% meshSize - mesh size (m)
% ilo - low x cell in PP line (cells)
% ihi - high x cell in PP line (cells)
% jlo - low y cell in PP line (cells)
% jhi - high y cell in PP line (cells)
% klo - low z cell in PP line (cells)
% khi - high z x cell in PP line (cells)
% fieldOrient - field orientation (1-6)
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
% Date: 28/08/2007

  processFileName = 'process.dat';

  fd = fopen( processFileName , 'rb' );
  
  if fd == -1
    error( 'error: cannot open file %s' , processFileName );
  end %if

  % Header.
  header = fgets( fd );
  
  % Number of OP lines.
  ko = fscanf( fd , '%d' , [ 1 ,1 ] );
  
  % OP lines.
  io = zeros( ko , 9 );
  for no=1:ko
      io(no,1:9) = fscanf( fd , '%d' , [ 1 , 9 ] );
  end % for
  
  % Start and stop frequency steps.
  m = fscanf( fd , '%d' , [ 1 , 2 ] );
  nt1 = m(1);
  nt2 = m(2);
  
  %if nt1 ~= nt2 
  %  warning( 'warning - only first frequency step will be output' );
  %end % if
  
  % Start, stop and frequency step.
  m = fscanf( fd , '%e' , [ 1 , 3 ] );
  d1 = m(1);
  d2 = m(2);
  ds = m(3);
  
  % Mesh size.
  meshSize = fscanf( fd , '%e' , [ 1 , 1 ] );
  
  % Processing points.
  m = fscanf( fd , '%d' , [ 1 , 7 ] );  
  ilo = m(1);
  ihi = m(2);
  jlo = m(3);
  jhi = m(4);
  klo = m(5);
  khi = m(6);
  fieldOrient = m(7);
  
  fclose( fd );

end % function
