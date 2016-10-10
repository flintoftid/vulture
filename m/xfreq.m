function xfreq( wfFileName , tdFileName , fdFileName )
%
% xfreq - FFT time domain ASCII observer and normalise to waveform.
%
% xfreq( wfFileName , tdFileName , fdFileName )
%
% Inputs:
%
% wfFileName  - name of time domain waveform file.
% tdFileName  - name of input time domain field observer file.
% fdFileName  - name of output frequency domain waveform file.
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
% Author: Ian Flintoft <ian.flintoft@googlemail.com>
%

  % Read reference waveform.
  data_wf = readDataFile( wfFileName );
 
  % Read time domain field data.
  data_td = readDataFile( tdFileName );

  % Sample times, time-step and number of time-steps.
  t = data_wf(:,2);
  dt = t(2) - t(1);
  N = length( t );

  % Waveform values.
  wf = data_wf(:,3);

  % Field values.
  ex = data_td(:,3);
  ey = data_td(:,4);
  ez = data_td(:,5);
  hx = data_td(:,6);
  hy = data_td(:,7);
  hz = data_td(:,8);

  % FFT and normalise to waveform.
  Wf = fft( wf );
  Ex = fft( ex ) ./ Wf;
  Ey = fft( ey ) ./ Wf;
  Ez = fft( ez ) ./ Wf;
  Hx = fft( hx ) ./ Wf;
  Hy = fft( hy ) ./ Wf;
  Hz = fft( hz ) ./ Wf;

  % FFT sample frequencies.
  f = (0:N-1).' ./ dt ./ N;

  % Reduce data.
  f = f(1:N/10);
  Ex = Ex(1:N/10);
  Ey = Ey(1:N/10);
  Ez = Ez(1:N/10);
  Hx = Hx(1:N/10);
  Hy = Hy(1:N/10);
  Hz = Hz(1:N/10);

  % Data table for output.
  data_fd = [ f , abs( Ex ) , angled( Ex ) , ...
                  abs( Ey ) , angled( Ey ) , ...
                  abs( Ez ) , angled( Ez ) , ...
                  abs( Hx ) , angled( Hx ) , ...
                  abs( Hy ) , angled( Hy ) , ...
                  abs( Hz ) , angled( Hz ) ];

  % Data file header.
  header = { 'f(MHz) ' , '|Ex| (V/m)' , '/_Ex (deg.)' , ...
                         '|Ey| (V/m)' , '/_Ey (deg.)' , ...
                         '|Ez| (V/m)' , '/_Ex (deg.)' , ...
                         '|Hx| (V/m)' , '/_Hx (deg.)' , ...
                         '|Hy| (V/m)' , '/_Hy (deg.)' , ...
                         '|Hz| (V/m)' , '/_Hz (deg.)' };

  % Write in double precision ASCII.
  writeDataFileDP( fdFileName , data_fd , header );
 
end %function
