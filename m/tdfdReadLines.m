function [ x , y , z ] = tdfdReadLines()
%
%  [ x , y , z ] = tdfdReadLines()
%
% Read in mesh lines from time domain code data files.
%
% Optimised for extracting data in 1D, 2D and 3D spatial slices at a small
% number of times. Not suitable for extracting long time
% slices at a single point.
% 
% Inputs:
%
% Outputs:
%
% x  - vector of x positions, x(i) (metres)
% y  - vector of y positions, y(j) (metres)
% z  - vector of z positions, z(k) (metres)
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
% Date: 04/11/2011

  % Names of the data files.
  xlinesFileName = 'xlines.dat';
  ylinesFileName = 'ylines.dat';
  zlinesFileName = 'zlines.dat';

  x = [];
  y = [];
  z = [];

  % Read in x lines.
  fdLines = fopen( xlinesFileName , 'r' );
  if fdLines == -1
    error( 'cannot open file %s' , xlinesFileName );
  end %if

  lineNum = 1;
  while( ~feof( fdLines ) )
    [ data , count ] = fscanf( fdLines , '%d %e' , 2 );
    if( count ~= 2 )
      if( feof( fdLines ) )
        break;
      else
        error( 'error reading from file %s' , xlinesFileName );
      end %if
    else
      ix(lineNum) = data( 1 );
      x(lineNum) = data( 2 );
      lineNum = lineNum + 1;
    end % if
  end % while

  fclose( fdLines );

  % Read in y lines.
  fdLines = fopen( ylinesFileName , 'r' );
  if fdLines == -1
    error( 'cannot open file %s' , ylinesFileName );
  end %if

  lineNum = 1;
  while( ~feof( fdLines ) )
    [ data , count ] = fscanf( fdLines , '%d %e' , 2 );
    if( count ~= 2 )
      if( feof( fdLines ) )
        break;
      else
        error( 'error reading from file %s' , ylinesFileName );
      end %if
    else
      jy(lineNum) = data( 1 );
      y(lineNum) = data( 2 );
      lineNum = lineNum + 1;
    end % if
  end % while

  fclose( fdLines );

  % Read in x lines.
  fdLines = fopen( zlinesFileName , 'r' );
  if fdLines == -1
    error( 'cannot open file %s' , zlinesFileName );
  end %if

  lineNum = 1;
  while( ~feof( fdLines ) )
    [ data , count ] = fscanf( fdLines , '%d %e' , 2 );
    if( count ~= 2 )
      if( feof( fdLines ) )
        break;
      else
        error( 'error reading from file %s' , zlinesFileName );
      end %if
    else
      kz(lineNum) = data( 1 );
      z(lineNum) = data( 2 );
      lineNum = lineNum + 1;
    end % if
  end % while

  fclose( fdLines );

end % function
  
