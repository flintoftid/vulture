function [ x , y , z , Ex , Ey , Ez , Hx , Hy , Hz , f ] = ...
    tdfdReadFrequencyDat0D( opCardNum , outputPhase , boundingBox )
%
%  [ x , y , z , Ex , Ey , Ez , Hx , Hy , Hz , f ] = ...
%    tdfdReadFrequencyDat0D( opCardNum , outputPhase , boundingBox )
%
% Read in frequency domain field magnitude and optionally
% phase data from time domain code binary data files.
%
% For extracting data at a single point over the whole frequency range.
% 
% Inputs:
%
% opCardNum   - which OP card number to extract data for.
%               Defaults to first card.
% outputPhase - boolean (0 or 1) determining whether or
%               not the phase data is read.
%               Default is not to read phase data.
% boundingBox - bounding box for position data to read in
%               form [ xlo , xhi , ylo , yhi, zlo , zhi ].
%               If not specified the values in the PP line
%               of the process.dat file are used.
%
% Outputs:
%
% x  - x position (cells)
% y  - y position (cells)
% z  - z position (cells)
%
% Ex - 1D array of x component of electric field, Ex(i,j,k,f) (V/m)
% Ey - 1D array of y component of electric field, Ey(i,j,k,f) (V/m)
% Ez - 1D array of z component of electric field, Ez(i,j,k,f) (V/m)
% Hx - 1D array of x component of magnetic field, Hx(i,j,k,f) (A/m)
% Hy - 1D array of y component of magnetic field, Hy(i,j,k,f) (A/m)
% Hz - 1D array of z component of magnetic field, Hz(i,j,k,f) (A/m)
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
% Date: 28/08/2007

  fprintf( '\n\n*** WARNING - Beta fast seek version - limited testing! ***\n\n' );

  % Define output arrays.
  x = [];
  y = [];
  z = [];
  Ex = [];
  Ey = [];
  Ez = [];
  Hx = [];
  Hy = [];
  Hz = [];

  % Process arguments.
  if nargin == 0
    opCardNum = 1;
    outputPhase = 0;
    boundingBox = [];
  elseif nargin == 1
    outputPhase = 0;
    boundingBox = [];
  elseif nargin == 2
    boundingBox = [];
  elseif nargin == 3
    ;
  else
    error( 'too many arguments' );
  end % if

  % Validate opCardNum.
  if opCardNum < 1
    error( 'invalid value for parameter opCardNum' );
  end % if

  % Validate phase request.
  if outputPhase ~= 0 && outputPhase ~= 1
      error( 'invalid value for parameter outputPhase' );
  end % if

  % Validate bounding box.
  if ~isempty( boundingBox )
    if length( boundingBox ) ~= 6 || any( boundingBox < 0 )
      error( 'invalid value for parameter boundingBox' );
    end % if
  end %if

  % Names of the data files.
  frequencyFileName = 'frequency.dat';
  phaseFileName = 'phase.dat';

  % Read the process.dat file.
  [ header , ko , io , nt1 , nt2 , d1 , d2 , ds , meshSize , ...
    ilo , ihi , jlo , jhi , klo , khi , fieldOrient ] =  tdfdReadProcessDat();

  % Bounding box of requested data
  if isempty( boundingBox )
    boundingBox = [ ilo , ihi , jlo , jhi , klo , khi ];
  end %if

  if ( boundingBox(1) ~=  boundingBox(2) ) || ...
     ( boundingBox(3) ~=  boundingBox(4) ) || ...
     ( boundingBox(6) ~=  boundingBox(6) )
    error( 'bounding box must be single point for 0D' );
  end % if

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
  freqNum=1:numFreqInFile;

  % Preallocate arrays.
  Ex = zeros( numFreqInFile , 1 );
  Ey = zeros( numFreqInFile , 1 );
  Ez = zeros( numFreqInFile , 1 );
  Hx = zeros( numFreqInFile , 1 );
  Hy = zeros( numFreqInFile , 1 );
  Hz = zeros( numFreqInFile , 1 );

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

  % Find offsets to required point in op block.
  opPreOffset = 0;
  found = 0;
  for kk=io(opCardNum,7):io(opCardNum,9):io(opCardNum,8)
    for jj=io(opCardNum,4):io(opCardNum,6):io(opCardNum,5)
      for ii=io(opCardNum,1):io(opCardNum,3):io(opCardNum,2)
        if ii == boundingBox(1) && jj == boundingBox(3) && kk == boundingBox(5)
          found = 1;
          break;
        end  %if
        opPreOffset = opPreOffset + 9 * 4;
      end % for
      if found 
        break;
      end %if
    end %for
    if found 
      break;
    end %if
  end %for
  if ~found
    error( 'point not found in OP block' );
  end %if
  opPostOffset = opBlockSize(opCardNum) - opPreOffset - 9 * 4;
  assert( opPreOffset + opPostOffset + 9 * 4 == opBlockSize(opCardNum) );

  % Make sure requested frequencies are in order
  freqNum = sort( freqNum );
  freq = ( freqNum - 1 ) * freqStep;
  freqSkip = [ freqNum(1) - 1 ,  diff( freqNum ) - 1 ];
  lastFreqSkip = numFreqInFile - freqNum(end);
  assert( sum( freqSkip ) + length( freqNum ) + lastFreqSkip == numFreqInFile );
  f = freqInFile( freqNum );

  % Read in data.
  for fq=1:length( freqNum )

    % Skip to required frequency block.
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

    % Skip to required point in op block.
    ret = fseek( fdFreq , opPreOffset , 'cof' );
    if ret ~= 0
      msg = ferror( fdFreq );
      error( msg );
    end %if
    if outputPhase == 1
      ret = fseek( fdPhase , opPreOffset , 'cof' );
      if ret ~= 0
        msg = ferror( fdPhase );
        error( msg );
      end %if
    end %if

    % Read data at point
    ix = fread( fdFreq , 1 , 'int32' );
    jy = fread( fdFreq , 1 , 'int32' );
    kz = fread( fdFreq , 1 , 'int32' );
    ex = fread( fdFreq , 1 , 'float32' );
    ey = fread( fdFreq , 1 , 'float32' );
    ez = fread( fdFreq , 1 , 'float32' );
    hx = fread( fdFreq , 1 , 'float32' );
    hy = fread( fdFreq , 1 , 'float32' );
    hz = fread( fdFreq , 1 , 'float32' );
    assert( ix == boundingBox(1) );
    assert( jy == boundingBox(3) );
    assert( kz == boundingBox(5) );

    if outputPhase == 1
      ix = fread( fdPhase , 1 , 'int32' );
      jy = fread( fdPhase , 1 , 'int32' );
      kz = fread( fdPhase , 1 , 'int32' );
      exp = fread( fdPhase , 1 , 'float32' );
      eyp = fread( fdPhase , 1 , 'float32' );
      ezp = fread( fdPhase , 1 , 'float32' );
      hxp = fread( fdPhase , 1 , 'float32' );
      hyp = fread( fdPhase , 1 , 'float32' );
      hzp = fread( fdPhase , 1 , 'float32' );
      assert( ix == boundingBox(1) );
      assert( jy == boundingBox(3) );
      assert( kz == boundingBox(5) );
    else
      exp = 0.0;
      eyp = 0.0;
      ezp = 0.0;
      hxp = 0.0;
      hyp = 0.0;
      hzp = 0.0;
    end % if outputPhase
    x = ix;
    y = jy;
    z = kz;
    Ex(fq) = ex * e^( j * exp / 180 * pi );
    Ey(fq) = ey * e^( j * eyp / 180 * pi );
    Ez(fq) = ez * e^( j * ezp / 180 * pi );
    Hx(fq) = hx * e^( j * hxp / 180 * pi );
    Hy(fq) = hy * e^( j * hyp / 180 * pi );
    Hz(fq) = hz * e^( j * hzp / 180 * pi );

    % Skip to end of op block.
    ret = fseek( fdFreq , opPostOffset , 'cof' );
    if ret ~= 0
      msg = ferror( fdFreq );
      error( msg );
    end %if
    if outputPhase == 1
      ret = fseek( fdPhase , opPostOffset , 'cof' );
      if ret ~= 0
        msg = ferror( fdPhase );
        error( msg );
      end %if
    end %if

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

  % Skip to end of data file(s) and check haved reached end of file.
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
