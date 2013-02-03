%============================================================================
%
% Constants class
% 
% This class contains general physical constants used in the calculations.
% For specific quantities, look in the other classes of the LIT package.
%
%============================================================================
classdef Constants
properties ( Access = public, Constant )          % with the 'Constant' attribute, they can be called from outside without an instance

    % conversion factors
    Na      = 6.02214179e23;                      % Avogadro number
    u_in_Kg = 1.660538782e-27;                    % unit of atomic mass
    kb      = 1.3806504e-23;                      % Boltzmann constant in J/K
    T0      = -273.15;                            % absolute zero in °C

    %end of properties
end

% end of LIT class
end
