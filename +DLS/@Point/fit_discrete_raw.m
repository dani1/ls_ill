function fit_obj = fit_discrete_raw ( t, g, dg, method, q, protein)
%============================================================================
% FIT CORRELOGRAMMS USING DISCRETE DECAYS
%============================================================================
% take correlogramm and Q-vector as input and return values of fit coefficient as output
% the fit is limited to discrete discrete times
    Min_D1   = 0.5;  % A^2 / ns % changed from 1 to 0.5
    Max_D1   = 100;  % A^2 / ns
    Start_D1 = 6;    % A^2 / ns
    Min_D2   = 0;    % A^2 / ns
    Max_D2   = 2.5;  % A^2 / ns
    Start_D2 = 6e-1; % A^2 / ns
    Gf       = @(x)1e6*x*q^2;    % convert limits for D into limits for Gammas
    Min_Gamma1   = Gf(Min_D1);   % ms^{-1}
    Max_Gamma1   = Gf(Max_D1);   % ms^{-1}
    Start_Gamma1 = Gf(Start_D1); % ms^{-1}
    Min_Gamma2   = Gf(Min_D2);   % ms^{-1}
    Max_Gamma2   = Gf(Max_D2);   % ms^{-1}
    Start_Gamma2 = Gf(Start_D2); % ms^{-1}
    %==  PART 1: set the fit parameters ==
    switch method % depending on the method chosen, se the parameters
    case 'SingleBeta'
        fit_function	= 'beta * (Ae * exp( - Gammae * t )) .^2';
        coeffnames = {'beta' 'Ae' 'Gammae'     };
        lowerbond  = [0      0    Min_Gamma2   ]; %changed 2 to Min_Gamma2 from Min_Gamma3
        upperbond  = [10     5    Max_Gamma1   ];
        startpoint = [1      1.0  Start_Gamma1 ];

    case 'DoubleFreeBeta' % double exponential decay
        fit_function	= 'beta *( A1 * exp( - Gamma1 * t ) + A2 * exp( - Gamma2 * t ) ).^2';
        coeffnames = {'beta' 'A1' 'Gamma1'     'A2' 'Gamma2'     };
        lowerbond  = [0      0    0            0    0            ];
        upperbond  = [10     1    Max_Gamma1   1    Max_Gamma1   ]; %
        startpoint = [1      0.1  Start_Gamma1 0.9  Start_Gamma2 ];
    otherwise
        error('Method not recognized!');
    end
    %==  PART 2: perform the fit ==
    switch method % depending on the method chosen, perform the fit
    case 'Cumulants'
        weights = 1 ./ dyc .^ 2;                                      % set the weights for the fit
        % set other fit options
        fit_options = fitoptions( 'Method', 'LinearLeastSquares', ...
        'Weights', weights );
        fit_type = fittype( fit_function,           ...
        'coefficients' , {'loga','gamma'},  ...
        'dependent'    , 'logg1',       ...
        'independent'  , 't',           ...
        'options'      , fit_options        );                        % A linear model is ok
        fit_obj = fit( t, yc, fit_type );
        otherwise
        weights = 1 ./ dg .^ 2;                                       % set the weights for the fit
        % set other fit options
        fit_options = fitoptions( 'Method','NonLinearLeastSquares',...
                    'Weights'     , weights          , ...
                    'MaxFunEvals' , 100000           , ...
                    'Lower'       , lowerbond        , ...
                    'Upper'       , upperbond        , ...
                    'Startpoint'  , startpoint              );

        fit_type = fittype(fit_function, 'independent','t', 	...
                    'options'      , fit_options     , ...
                    'coefficients' , coeffnames         );
        % perform the numerical fit
        % TODO: control the behaviour of the correlations at 100 times higher than the
        %        first decay time. If there is still something important, do not fit
        %        the data and impose NaN for both Gamma and dGamma. One can not live with
        %        scorpions is his bed!
        fit_obj = fit( t, g, fit_type );
    end
end	% fit_discrete

