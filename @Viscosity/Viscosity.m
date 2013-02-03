%============================================================================
%
% Viscosity class
% 
% This class contains the viscosity of various substances, used in the data
% analysis
% 
% vim: foldmethod=marker
%
%============================================================================
classdef Viscosity

 % Properties {{{1

properties ( Access = public, Constant )    % with the 'Constant' attribute, they can be called from outside without an instance

    water_Units = 'cP';

%end of properties
end

% Methods {{{1
methods ( Static )

  % water viscosity as a function of temperature, in centiPoise {{{2
    function eta = water( T )

        T0 = 225.334;   % K
        A  = 802.25336;  % cP * K^g
        a  = 3.4741e-3;  % K^-1
        b  = -1.7413e-5; % K^-2
        c  = 2.7719e-8;  % K^-3
        g  = 1.53026;
    
        if T < 100
            T = T - LIT.Constants.T0;
        end
        DT  = T -T0; % Delta T
        eta = A .* ( DT + a*DT.^2 + b*DT.^3 + c*DT.^4 ).^(-g);
    end % }}}2

end

% end of LIT class
end
