classdef Point
% This class describes a single SLS datapoint

properties

    Protein % what protein?
    Salt % what salt?
    C % protein concentration
    Unit_C   = 'g/l'
    Cs % salt concentration
    Unit_Cs  = 'mM'

    Angle
    Unit_Q   = 'A^{-1}'

    Unit_KcR = 'mol g^{-1}'
    Unit_X_T = 'l * J^{-1}';

    T
    Unit_T   = 'K'
    n
    dndc

end

properties ( Dependent, SetAccess = private )

    Phi                       % volume fraction
    Q
    KcR
    dKcR
    X_T
    dX_T

end

properties ( Hidden )

    Instrument
    C_set
    n_set
    dndc_set
    KcR_raw
    dKcR_raw
    datetime
    datetime_raw
    Imon
    dImon

end

methods

    function Phi = get.Phi ( self )
        v0  = LIT.(self.Protein).v0;
        Phi = self.C * v0;
    end

    function Q = get.Q ( self )
        Qexp = '4 * pi * n * sind( 0.5 * theta ) / lambda';
        Qf   = inline(Qexp,'n','theta','lambda');
        Q    = Qf( self.n, self.Angle, self.Instrument.Lambda );
    end

    function KcR = get.KcR ( self )
        KcRexpr	= 'old * ( c / c_set ) * ( dndc / dndc_set )^2 * ( n / n_set )^2';
        KcRf = inline(KcRexpr,'old',        ...
                    'c'   , 'c_set'       , ...
                    'dndc', 'dndc_set'    , ...
                    'n'   , 'n_set'     );
        KcR  = KcRf( self.KcR_raw,           ...
                    self.C   , self.C_set    , ...
                    self.dndc, self.dndc_set , ...
                    self.n   , self.n_set  );
    end

    function dKcR = get.dKcR ( self )
        dKcR = self.KcR / self.KcR_raw * self.dKcR_raw;
    end

    function X_T = get.X_T ( self )
    % the osmotic isothermal compressibility is related to the structure factor at
    % vanishing Q by the following formula:
    %
    % X :=   1/c * ( dc / dP )_T = 1 / [ Na * kT * c * M * S(c, q -> 0)  ]
                                %= 1 / [ Na * kT * c *   (Kc/R)      ]
    %
    % where X is the compressibility, c is the concentration in mass/volume, N_A is
    % Avogadro's number and S(c,q) is the static structure factor.
    %
    % Please note that X(c) is independent on Q (because it is defined for Q-->0) and on M.

        X_T	= SLS.X_Tf( self.T, self.C, self.KcR );
    end

    function dX_T = get.dX_T ( self )
        dX_T = self.X_T .* self.dKcR ./ self.KcR;            % TODO: errors only come from KcR
    end

end

end
