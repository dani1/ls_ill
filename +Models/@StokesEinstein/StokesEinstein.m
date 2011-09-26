%============================================================================
%
% Stokes-Einstein class
% 
% This class contains a couple of functions to use for the SE relation
% 
% vim: foldmethod=marker
%
%============================================================================
classdef StokesEinstein

 % Properties {{{1
 properties ( Access = public, Constant )	% with the 'Constant' attribute, they can be called from outside without an instance

  formula	= 'D0 = k*T / ( 6 * pi * eta0 * rh )'
  dilutelimit_diffusionconstant_Units = 'A^2 / ns'
  hydrodynamic_radius_Units = 'm'

 end % }}}1
 % Methods {{{1
 methods ( Static )

  % D0 as a function of the hydrodynamic radius, eta, and T {{{2
  function D0 = dilutelimit_diffusionconstant( rh, eta0, T )

   k = LIT.Constants.kb;
   if T < 100
    T = T - LIT.Constants.T0;
   end

   D0	= 1e14 * k * T ./ ( 6 * pi * eta0 * rh );

  end % }}}2
  % Hydrodynamic radius as a function of D0, eta, and T {{{2
  function rh = hydrodynamic_radius( D0, eta0, T )

   k = LIT.Constants.kb;
   if T < 100
    T = T - LIT.Constants.T0;
   end

   rh	= 1e14 * k * T ./ ( 6 * pi * eta0 * D0 );

  end % }}}2

 end % }}}1

% end of LIT class
end
