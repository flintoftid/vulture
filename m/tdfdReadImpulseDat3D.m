function [ x , y , z , Ex , Ey , Ez , Hx , Hy , Hz , t ] = ...
    tdfdReadImpulseDat3D( opCardNum , timeNum )
%
%  [ x , y , z , Ex , Ey , Ez , Hx , Hy , Hz , t ] = ...
%    tdfdReadImpulseDat3D( opCardNum ,  timeNum )
%
% Read in time domain field data from time domain code binary data files.
%
% Optimised for extracting data full data-set or 1D, 2D and 3D spatial 
% slices at a small number of times.
%
% Inputs:
%
% opCardNum   - which OP card number to extract data for.
%               Defaults to first card.
% timeNum     - vector of time step numbers to extract
%               If not specified the range nt1:nt2 in the
%               process.dat file are used.
%
% Outputs:
%
% x  - vector of x positions, x(i) (cells)
% y  - vector of y positions, y(j) (cells)
% z  - vector of z positions, z(k) (cells)
%
% Ex - 4D array of x component of electric field, Ex(i,j,k,t) (V/m)
% Ey - 4D array of y component of electric field, Ey(i,j,k,t) (V/m)
% Ez - 4D array of z component of electric field, Ez(i,j,k,t) (V/m)
% Hx - 4D array of x component of magnetic field, Hx(i,j,k,t) (A/m)
% Hy - 4D array of y component of magnetic field, Hy(i,j,k,t) (A/m)
% Hz - 4D array of z component of magnetic field, Hz(i,j,k,t) (A/m)
%
% t  - vector of times of extracted data (s)
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

  % Process arguments.
  if nargin == 0
    opCardNum = 1;
    timeNum = [];
    boundingBox = [];
  elseif nargin == 1
    timeNum = [];
    boundingBox = [];
  elseif nargin == 2
    ;
  else
    error( 'too many arguments' );
  end % if

  % Validate requested frequency number vector
  if ~isempty( timeNum )
    if any( timeNum ) < 0 
      error( 'invalid value for parameter timeNum' );
    end % if
  end %if

  % Validate opCardNum.
  if opCardNum < 1
    error( 'invalid value for parameter opCardNum' );
  end % if

  % Names of the data files.
  impulseFileName = 'impulse.dat';

  % Read the process.dat file.
  [ header , ko , io , nt1 , nt2 , d1 , d2 , ds , meshSize , ...
    ilo , ihi , jlo , jhi , klo , khi , fieldOrient ] =  tdfdReadProcessDat();

  % Check requested OP card exists.
  if opCardNum > ko
    error( 'opCardNum greater than number of OP cards in data set' );
  end % if

  % Open data files.
  fdImp = fopen( impulseFileName , 'rb' );
  if fdImp == -1
    error( 'cannot open file %s' , impulseFileName );
  end %if

  % Read in number of time steps and mesh size.
  numTimeInFile = fread( fdImp , 1 , 'int32' );
  timeStep = fread( fdImp , 1 , 'float32' );
  timeStepInFile = timeStep .* (0:(numTimeInFile-1));

  fprintf( 'Number of time-steps in file = %d.\n' , numTimeInFile );
  fprintf( 'Time step = %g s\n' , timeStep );

  % If frequency numbers not specified choose all.
  if isempty( timeNum )
    timeNum=1:numTimeInFile;
  elseif( any( timeNum > numTimeInFile ) )
    error( 'requested time step number too large for data file!' );
  end % if

  % Determine size of each OP block and frequency block.
  for no=1:ko
    opBlockSize(no) = 9 * 4 * length( io(no,7):io(no,9):io(no,8) ) ...
                            * length( io(no,4):io(no,6):io(no,5) ) ...
                            * length( io(no,1):io(no,3):io(no,2) );
  end % for

  timeBlockSize = sum( opBlockSize );
  totalFileSize = 8 + numTimeInFile * timeBlockSize;
  opPreSkip = sum( opBlockSize(1:opCardNum-1) );
  opPostSkip = sum( opBlockSize(opCardNum+1:end) );

  assert( opPreSkip + opPostSkip + opBlockSize(opCardNum) == timeBlockSize );

  % Make sure requested times are in order
  timeNum = sort( timeNum );
  time = ( timeNum - 1 ) * timeStep;
  timeSkip = [ timeNum(1) - 1 ,  diff( timeNum ) - 1 ];
  lastTimeSkip = numTimeInFile - timeNum(end);
  assert( sum( timeSkip ) + length( timeNum ) + lastTimeSkip == numTimeInFile );
  t = timeStepInFile( timeNum );

  % Number of output points in each direction.
  nx = length( io(opCardNum,1):io(opCardNum,3):io(opCardNum,2) );
  ny = length( io(opCardNum,4):io(opCardNum,6):io(opCardNum,5) );
  nz = length( io(opCardNum,7):io(opCardNum,9):io(opCardNum,8) );

  % Initialise arrays.
  Ex = zeros( nx , ny , nz , length( timeNum ) );
  Ey = zeros( nx , ny , nz , length( timeNum ) );
  Ez = zeros( nx , ny , nz , length( timeNum ) );
  Hx = zeros( nx , ny , nz , length( timeNum ) );
  Hy = zeros( nx , ny , nz , length( timeNum ) );
  Hz = zeros( nx , ny , nz , length( timeNum ) );
  x = io(opCardNum,1):io(opCardNum,3):io(opCardNum,2);
  y = io(opCardNum,4):io(opCardNum,6):io(opCardNum,5);
  z = io(opCardNum,7):io(opCardNum,9):io(opCardNum,8);

  % Read in data.
  for ts=1:length( timeNum )

    % Skip to required time block.
    % fprintf( 'Skipping to time step number %d (%g ns)...\n' , timeNum(ts) , time(ts) / 1e-9 );

    if( timeSkip( ts ) ~= 0  )

      ret = fseek( fdImp , timeSkip( ts ) * timeBlockSize , 'cof' );
      if ret ~= 0
        msg = ferror( fdImp );
        error( msg );
      end %if

    end %if 

    % Skip to required OP block.
    ret = fseek( fdImp , opPreSkip , 'cof' );
    if ret ~= 0
      msg = ferror( fdImp );
      error( msg );
    end %if

    % Read required OP block.
    % fprintf( 'Reading data for time step number %d (%g ns)...\n' , timeNum(ts) , time(ts) / 1e-9 );

    % Save start of OP block.
    pos = ftell( fdImp );

    % Read entire block as vector of integers and reshape into array.
    tmpInd = fread( fdImp , 9 * nx * ny * nz , 'int32' ); 
    tmpInd2 = reshape( tmpInd , [ 9 , nx , ny , nz ] );
    tmpInd = permute( tmpInd2 , [ 2 , 3 , 4 , 1 ] );

    % Extract required indices and verify.
    xx = tmpInd(:,1,1,1);
    yy = tmpInd(1,:,1,2);
    zz = tmpInd(1,1,:,3);

    assert( xx(:) == x(:) );
    assert( yy(:) == y(:) );
    assert( zz(:) == z(:) );

    % Return to start of OP block.
    ret = fseek( fdImp , pos );
    if ret ~= 0
      msg = ferror( fdImp );
      error( msg );
    end %if

    % Read entire block as vector of floats and reshape into array.
    tmpField = fread( fdImp , 9 * nx * ny * nz , 'float32' ); 
    tmpField2 = reshape( tmpField , [ 9 , nx , ny , nz ] );
    tmpField = permute( tmpField2 , [ 2 , 3 , 4 , 1 ] );

    % Extract fields into arrays.
    Ex(:,:,:,ts) = tmpField(:,:,:,4);
    Ey(:,:,:,ts) = tmpField(:,:,:,5);
    Ez(:,:,:,ts) = tmpField(:,:,:,6);
    Hx(:,:,:,ts) = tmpField(:,:,:,7);
    Hy(:,:,:,ts) = tmpField(:,:,:,8);
    Hz(:,:,:,ts) = tmpField(:,:,:,9);

    % fprintf( 'Extracted data for time step number %d (%g ns)...\n' , timeNum(ts) , time(ts) / 1e-9 );

    % Skip remaining OP blocks in file(s).
    ret = fseek( fdImp , opPostSkip , 'cof' );
    if ret ~= 0
      msg = ferror( fdImp );
      error( msg );
    end %if

  end %for ts

  % Skip to end of data file and check haved reached end of file.
  % fprintf( 'Skipping to end of file...\n' );

  ret = fseek( fdImp , lastTimeSkip * timeBlockSize , 'cof' );
  if ret ~= 0
    msg = ferror( fdImp );
    error( msg );
  end %if

  assert( ftell( fdImp ) == totalFileSize );

  % Close data file(s).
  fclose( fdImp );

end % function

