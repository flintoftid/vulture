function [ x , y , z , Ex , Ey , Ez , Hx , Hy , Hz , f ] = ...
    tdfdReadFrequencyDat3D( opCardNum , freqNum , outputPhase )
%
%  [ x , y , z , Ex , Ey , Ez , Hx , Hy , Hz , f ] = ...
%    tdfdReadFrequencyDat3D( opCardNum , freqNum , outputPhase )
%
% Read in frequency domain field magnitude and optionally
% phase data from time domain code binary data files.
%
% Optimised for extracting fuill data-set and 1D, 2D and 3D spatial slices 
% at a small number of frequency points. 
% 
% Inputs:
%
% opCardNum   - which OP card number to extract data for.
%               Defaults to first card.
% freqNum     - vector of frequency numbers to extract
%               If not specified the range nt1:nt2 in the
%               process.dat file are used.
% outputPhase - boolean (0 or 1) determining whether or
%               not the phase data is read.
%               Default is not to read phase data.
%
% Outputs:
%
% x  - vector of x positions, x(i) (cells)
% y  - vector of y positions, y(j) (cells)
% z  - vector of z positions, z(k) (cells)
%
% Ex - 4D array of x component of electric field, Ex(i,j,k,f) (V/m)
% Ey - 4D array of y component of electric field, Ey(i,j,k,f) (V/m)
% Ez - 4D array of z component of electric field, Ez(i,j,k,f) (V/m)
% Hx - 4D array of x component of magnetic field, Hx(i,j,k,f) (A/m)
% Hy - 4D array of y component of magnetic field, Hy(i,j,k,f) (A/m)
% Hz - 4D array of z component of magnetic field, Hz(i,j,k,f) (A/m)
%
% f  - vector of frequencies of extracted data (Hz)
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
    freqNum = [];
    outputPhase = 0;
  elseif nargin == 1
    freqNum = [];
    outputPhase = 0;
  elseif nargin == 2
    outputPhase = 0;
  elseif nargin == 3
    ;
  else
    error( 'too many arguments' );
  end % if

  % Validate requested frequency number vector
  if ~isempty( freqNum )
    if any( freqNum ) < 0 
      error( 'invalid value for parameter freqNum' );
    end % if
  end %if

  % Validate opCardNum.
  if opCardNum < 1
    error( 'invalid value for parameter opCardNum' );
  end % if

  % Validate phase request.
  if outputPhase ~= 0 && outputPhase ~= 1
      error( 'invalid value for parameter outputPhase' );
  end % if

  % Names of the data files.
  frequencyFileName = 'frequency.dat';
  phaseFileName = 'phase.dat';

  % Read the process.dat file.
  [ header , ko , io , nt1 , nt2 , d1 , d2 , ds , meshSize , ...
    ilo , ihi , jlo , jhi , klo , khi , fieldOrient ] =  tdfdReadProcessDat();

  % Check requested OP card exists.
  if opCardNum > ko
    error( 'opCardNum greater than number of OP cards in data set' );
  end % if

  % Open data files.
  fdFreq = fopen( frequencyFileName , 'rb' );
  if fdFreq == -1
    error( 'cannot open file %s' , frequencyFileName );
  end %if

  if outputPhase == 1
    fdPhase = fopen(phaseFileName , 'rb' );
    if fdPhase == -1
      error( 'cannot open file %s' , phaseFileName );
    end %if
  end % if

  % Read in number of frequencies and frequency step.
  numFreqInFile = fread( fdFreq , 1 , 'int32' );
  freqStep = fread( fdFreq , 1 , 'float32' );
  freqInFile = freqStep .* (0:(numFreqInFile-1));

  if outputPhase == 1
    numFreqInFile2 = fread( fdPhase , 1 , 'int32' );
    freqStep2 = fread( fdPhase , 1 , 'float32' );
    if ( numFreqInFile2 ~= numFreqInFile ) || ( freqStep2 ~= freqStep )
      error( 'size mismatch between frequency and phase data files' );
    end %if
  end % if

  fprintf( 'Number of frequencies in file = %d.\n' , numFreqInFile );
  fprintf( 'Frequency step = %g MHz.\n' , freqStep / 1e6 );

  % If frequency numbers not specified choose all.
  if isempty( freqNum )
    freqNum=1:numFreqInFile;
  elseif( any( freqNum > numFreqInFile ) )
    error( 'requested frequency number too large for data file!' );
  end % if

  % Determine size of each OP block and frequency block.
  for no=1:ko
    opBlockSize(no) = 9 * 4 * length( io(no,7):io(no,9):io(no,8) ) ...
                            * length( io(no,4):io(no,6):io(no,5) ) ...
                            * length( io(no,1):io(no,3):io(no,2) );
  end % no

  freqBlockSize = sum( opBlockSize );
  totalFileSize = 8 + numFreqInFile * freqBlockSize;
  opPreSkip = sum( opBlockSize(1:opCardNum-1) );
  opPostSkip = sum( opBlockSize(opCardNum+1:end) );

  assert( opPreSkip + opPostSkip + opBlockSize(opCardNum) == freqBlockSize );

  %fprintf( 'Frequency block size = %d.\n' , freqBlockSize );
  %fprintf( 'Total file size = %d.\n' , totalFileSize );
  %fprintf( 'OP pre-skip = %d.\n' , opPreSkip );
  %fprintf( 'OP post-skip = %d.\n' , opPostSkip );

  % Make sure requested frequencies are in order
  freqNum = sort( freqNum );
  freq = ( freqNum - 1 ) * freqStep;
  freqSkip = [ freqNum(1) - 1 ,  diff( freqNum ) - 1 ];
  lastFreqSkip = numFreqInFile - freqNum(end);
  assert( sum( freqSkip ) + length( freqNum ) + lastFreqSkip == numFreqInFile );
  f = freqInFile( freqNum );

  % Number of output points in each direction.
  nx = length( io(opCardNum,1):io(opCardNum,3):io(opCardNum,2) );
  ny = length( io(opCardNum,4):io(opCardNum,6):io(opCardNum,5) );
  nz = length( io(opCardNum,7):io(opCardNum,9):io(opCardNum,8) );
  %fprintf( 'OP shape = [%d,%d,%d].\n' , nx , ny , nz );

  % Initialise arrays.
  Ex = zeros( nx , ny , nz , length( freqNum ) );
  Ey = zeros( nx , ny , nz , length( freqNum ) );
  Ez = zeros( nx , ny , nz , length( freqNum ) );
  Hx = zeros( nx , ny , nz , length( freqNum ) );
  Hy = zeros( nx , ny , nz , length( freqNum ) );
  Hz = zeros( nx , ny , nz , length( freqNum ) );
  x = io(opCardNum,1):io(opCardNum,3):io(opCardNum,2);
  y = io(opCardNum,4):io(opCardNum,6):io(opCardNum,5);
  z = io(opCardNum,7):io(opCardNum,9):io(opCardNum,8);

  % Read in data.
  for fq=1:length( freqNum )

    % Skip to required frequency block.
    % fprintf( 'Skipping to frequency number %d (%g MHz)...\n' , freqNum(fq) , freq(fq) / 1e6 );

    if( freqSkip( fq ) ~= 0  )

      ret = fseek( fdFreq , freqSkip( fq ) * freqBlockSize , 'cof' );
      if ret ~= 0
        msg = ferror( fdFreq );
        error( msg );
      end %if

      if outputPhase == 1
        ret = fseek( fdPhase , freqSkip( fq ) * freqBlockSize , 'cof' );
        if ret ~= 0
          msg = ferror( fdPhase );
          error( msg );
        end %if
      end % if

    end %if 

    % Skip to required OP block.
    ret = fseek( fdFreq , opPreSkip , 'cof' );
    if ret ~= 0
      msg = ferror( fdFreq );
      error( msg );
    end %if

    if outputPhase == 1
      ret = fseek( fdPhase , opPreSkip , 'cof' );
      if ret ~= 0
        msg = ferror( fdPhase );
        error( msg );
      end %if
    end % if

    % Read required OP block.
    % fprintf( 'Reading data for frequency number %d (%g MHz)...\n' , freqNum(fq) , freq(fq) / 1e6 );

    % Save start of OP block.
    freqPos = ftell( fdFreq );

    % Read entire block as vector of integers and reshape into array.
    tmpInd = fread( fdFreq , 9 * nx * ny * nz , 'int32' ); 
    tmpInd2 = reshape( tmpInd , [ 9 , nx , ny , nz ] );
    tmpInd = permute( tmpInd2 , [ 2 , 3 , 4 , 1 ] );

    % Extract required indices and verify.
    xx = tmpInd(:,1,1,1);
    yy = tmpInd(1,:,1,2);
    zz = tmpInd(1,1,:,3);

    assert( all( xx(:) == x(:) ) );
    assert( all( yy(:) == y(:) ) );
    assert( all( zz(:) == z(:) ) );

    % Return to start of OP block.
    ret = fseek( fdFreq , freqPos , 'bof' );
    if ret ~= 0
      msg = ferror( fdFreq );
      error( msg );
    end %if

    % Read entire block as vector of floats and reshape into array.
    tmpAmp = fread( fdFreq , 9 * nx * ny * nz , 'float32' ); 
    tmpAmp2 = reshape( tmpAmp , [ 9 , nx , ny , nz ] );
    tmpAmp = permute( tmpAmp2 , [ 2 , 3 , 4 , 1 ] );

    if outputPhase == 1

      % Save start of OP block.
      phasePos = ftell( fdPhase );

      % Read entire block as vector of integers and reshape into array.
      tmpInd = fread( fdPhase , 9 * nx * ny * nz , 'int32' ); 
      tmpInd2 = reshape( tmpInd , [ 9 , nx , ny , nz ] );
      tmpInd = permute( tmpInd2 , [ 2 , 3 , 4 , 1 ] );

      % Extract required indices and verify.
      xx = tmpInd(:,1,1,1);
      yy = tmpInd(1,:,1,2);
      zz = tmpInd(1,1,:,3);

      assert( all( xx(:) == x(:) ) );
      assert( all( yy(:) == y(:) ) );
      assert( all( zz(:) == z(:) ) );

      % Return to start of OP block.
      ret = fseek( fdPhase , phasePos , 'bof' );
      if ret ~= 0
        msg = ferror( fdPhase );
        error( msg );
      end %if

      % Read entire block as vector of floats and reshape into array.
      tmpPhase = fread( fdPhase , 9 * nx * ny * nz , 'float32' ); 
      tmpPhase2 = reshape( tmpPhase , [ 9 , nx , ny , nz ] );
      tmpPhase = permute( tmpPhase2 , [ 2 , 3 , 4 , 1 ] );

      Ex(:,:,:,fq) = tmpAmp(:,:,:,4) .* exp( j .* tmpPhase(:,:,:,4) ./ 180 .* pi );
      Ey(:,:,:,fq) = tmpAmp(:,:,:,5) .* exp( j .* tmpPhase(:,:,:,5) ./ 180 .* pi );
      Ez(:,:,:,fq) = tmpAmp(:,:,:,6) .* exp( j .* tmpPhase(:,:,:,6) ./ 180 .* pi );
      Hx(:,:,:,fq) = tmpAmp(:,:,:,7) .* exp( j .* tmpPhase(:,:,:,7) ./ 180 .* pi );
      Hy(:,:,:,fq) = tmpAmp(:,:,:,8) .* exp( j .* tmpPhase(:,:,:,8) ./ 180 .* pi );
      Hz(:,:,:,fq) = tmpAmp(:,:,:,9) .* exp( j .* tmpPhase(:,:,:,9) ./ 180 .* pi );

    else

      Ex(:,:,:,fq) = tmpAmp(:,:,:,4);
      Ey(:,:,:,fq) = tmpAmp(:,:,:,5);
      Ez(:,:,:,fq) = tmpAmp(:,:,:,6);
      Hx(:,:,:,fq) = tmpAmp(:,:,:,7);
      Hy(:,:,:,fq) = tmpAmp(:,:,:,8);
      Hz(:,:,:,fq) = tmpAmp(:,:,:,9);

    end % if

    % fprintf( 'Extracted data for frequency number %d (%g MHz)...\n' , freqNum(fq) , freq(fq) / 1e6 );

    % Skip remaining OP blocks in file(s).
    ret = fseek( fdFreq , opPostSkip , 'cof' );
    if ret ~= 0
      msg = ferror( fdFreq );
      error( msg );
    end %if

    if outputPhase == 1
      ret = fseek( fdPhase , opPostSkip , 'cof' );
      if ret ~= 0
        msg = ferror( fdPhase );
        error( msg );
      end %if
    end % if

  end %for fq

  % Skip to end of data file(s) and check have reached end of file.
  % fprintf( 'Skipping to end of file...\n' );

  ret = fseek( fdFreq , lastFreqSkip * freqBlockSize , 'cof' );
  if ret ~= 0
    msg = ferror( fdFreq );
    error( msg );
  end %if
  assert( ftell( fdFreq ) == totalFileSize );

  if outputPhase == 1
    ret = fseek( fdPhase , lastFreqSkip * freqBlockSize , 'cof' );
    if ret ~= 0
      msg = ferror( fdPhase );
      error( msg );
    end %if
    assert( ftell( fdPhase ) == totalFileSize );
  end %if

  % Close data file(s).
  fclose( fdFreq );
  if outputPhase == 1
    fclose( fdPhase );
  end  % if

end % function
