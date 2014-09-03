function [ dBValue ] = db20( linearValue , minDBvalue )
%
% Convert linear field quantity to decibels.
%

  if nargin == 1
    minDBvalue = -500;
  end % if

  minLinearValue = ones( size(linearValue) ) * 10^( minDBvalue / 20.0 );

  dBValue = 20.0 * log10( abs( max( minLinearValue , linearValue ) ) );

end % function

