function tdfdWriteProcessDat( header , ko , io , nt1 , nt2 , d1 , d2 , ...
                              ds , meshSize , ilo , ihi , jlo , ...
                              jhi , klo , khi , fieldOrient )
%
% tdfdWriteProcessDat( ko , io , nt1 , nt2 , d1 , d2 , ds , ...
%                      meshSize , ilo , ihi , jlo , jhi , ...
%                      klo , khi , fieldOrient )
%
% Write processing control data from time domain code to
% process.dat file.
%
% Inputs:
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
% Outputs:
% 
% None.
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

  processFileName = 'process.dat';

  % Backup existing process.dat file
  copyfile( processFileName , [ processFileName , '~' ] );
  
  fd = fopen( processFileName , 'wb' );
  
  if fd == -1
    error( 'error: cannot open file %s' , processFileName );
  end %if

  % Header.
  fprintf( fd , '%s' , header );
  
  % Number of OP lines.
  fprintf( fd , '%d\n' , ko );
  
  % OP lines.
  for no=1:ko
      fprintf( fd , '%d ' , io(no,:) );
      fprintf( fd , '\n' );
  end % for

  % Start and stop frequency steps.
  fprintf( fd , '%d %d\n' , nt1 , nt2 );

  % Start, stop and frequency step.
  fprintf( fd , '%g %g %g\n' , d1 , d2 , ds );
  
  % Mesh size.
  fprintf( fd , '%g\n' , meshSize );
  
  % Processing points.
  fprintf( fd , '%d %d %d %d %d %d %d\n' , ...
           ilo , ihi , jlo , jhi , klo , khi , ...
           fieldOrient );
  
  fclose( fd );

end % function
