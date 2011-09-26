%============================================================================
% FIT CORRELOGRAMMS USING DISCRETE DECAYS
%============================================================================
function [ coeffnames coeffval dcoeffval ] = fit_discrete ( t, g, dg, method, q, system)
% take correlogramm and Q-vector as input and return values of fit coefficient as output
% the fit is limited to discrete discrete times

  Min_D1		= 1;			% A^2 / ns
  Max_D1		= 100;			% A^2 / ns
  Start_D1		= LIT.(system).D0;	% A^2 / ns
  Min_D2		= 0;			% A^2 / ns
  Max_D2		= 2.5;			% A^2 / ns
  Start_D2		= 10*LIT.(system).D0;	% A^2 / ns

  Min_Gamma1		= 1e6*Min_D1*q^2;	% ms^{-1}
  Max_Gamma1		= 1e6*Max_D1*q^2;	% ms^{-1}
  Start_Gamma1		= 1e6*Start_D1*q^2;	% ms^{-1}
  Min_Gamma2		= 1e6*Min_D2*q^2;	% ms^{-1}
  Max_Gamma2		= 1e6*Max_D2*q^2;	% ms^{-1}
  Start_Gamma2		= 1e6*Start_D2*q^2;	% ms^{-1}

 %==  PART 1: set the fit parameters ==

 switch method									% depending on the method chosen, se the parameters

  case 'Single'												% single exponential decay
   fit_function	= 'Ae * exp( - 2 * Gammae * t )';
   coeffnames	= { 'Ae'	'Gammae'	};
   lowerbond	= [ 0.8		Min_Gamma1	];
   upperbond	= [ 1.2		Max_Gamma1	];
   startpoint	= [ 1.0		Start_Gamma1	];

   case 'Streched'											% single streched decay
   fit_function	= 'As * exp( - 2 * (Gammas * t)^b )';
   coeffnames	= { 'As'	'Gammas'	'b'	};
   lowerbond	= [ 0.8		Min_Gamma1 	0	];
   upperbond	= [ 1.2		Max_Gamma1 	1	];
   startpoint	= [ 1.1		Start_Gamma1	0.8	];

  case 'Double'												% double exponential decay
   fit_function	= '( A1 * exp( - Gamma1 * t ) + A2 * exp( - Gamma2 * t ) ).^2';
   coeffnames	= { 'A1'	'Gamma1'	'A2'	'Gamma2'	};
   lowerbond	= [ 0.5		Min_Gamma1 	0	Min_Gamma2 	];
   upperbond	= [ 1		Max_Gamma1 	0.5	Max_Gamma2 	];
   startpoint	= [ 0.9		Start_Gamma1	0.1	Start_Gamma2	];

  case 'Single+Streched'										% double decay, exp + streched
   fit_function	= '( Ae * exp( - Gammae * t ) + As * exp(-( Gammas *t).^b )) .^2';
   coeffnames	= { 'Ae'	'Gammae'	'As'	'Gammas'	'b'	};
   lowerbond	= [ 0.5		Min_Gamma1 	0	Min_Gamma2 	0	];
   upperbond	= [ 1		Max_Gamma1 	0.5	Max_Gamma2 	1	];
   startpoint	= [ 0.9		Start_Gamma1	0.1	Start_Gamma2	0.8	];

  case 'DoubleStreched'										% double decay, exp + streched
   fit_function	= '( As1 * exp( - ( Gammas1 * t).^b1 ) + As2 * exp(-( Gammas2 *t).^b2 )) .^2';
   coeffnames	= { 'As1'	'Gammas1'	'As2'	'Gammas2'	'b1'	'b2'	};
   lowerbond	= [ 0.5		Min_Gamma1 	0	Min_Gamma2 	0	0	];
   upperbond	= [ 1		Max_Gamma1 	0.5	Max_Gamma2 	1	1	];
   startpoint	= [ 0.9		Start_Gamma1	0.1	Start_Gamma2	0.8	0.8	];

  otherwise
   error('Method not recognized!');
 end

 %==  PART 2: perform the fit ==

 weights =  1 ./ dg .^ 2;					% set the weights for the fit

 % set other fit options
 fit_options = fitoptions( 'Method','NonLinearLeastSquares',...
    		'Weights', weights,			...
    		'MaxFunEvals',100000, 			...
    		'Lower',lowerbond, 			...
    		'Upper',upperbond, 			...
    		'Startpoint',startpoint				);
 
 fit_type = fittype(fit_function, 'independent','t', 	...
    		'options', fit_options,			...
    		'coefficients',	coeffnames			);

 % perform the numerical fit
 % TODO: control the behaviour of the correlations at 100 times higher than the
 %        first decay time. If there is still something important, do not fit
 %        the data and impose NaN for both Gamma and dGamma. One can not live with
 %        scorpions is his bed!
 [ cf gof ] = fit( t, g, fit_type );

 % get fit result
 confidence	= confint(cf,0.95);
 coeffval	= coeffvalues (cf); 
 dcoeffval	=  0.5*(confidence(2,:)-confidence(1,:));

 %==  PART 3: output coefficients ==

 % add also diffusion constants for gamma coeffnames (decay rates)
 for i = 1 : length(coeffnames)
  if ~isempty(regexp(coeffnames{i},'Gamma'))		% if the parameter is a decay...

   Dname	= regexprep(coeffnames{i},'Gamma','D');	% ...set the name of the diffusion constant...
   Dval		= 1e-6 * coeffval(i) / q^2;		% ...set the diffusion constant...
   dDval	= 1e-6 * dcoeffval(i) / q^2;		% ...and its error...

   coeffval	= [ coeffval		Dval ];		% ...add D to the coeffval vector...
   dcoeffval	= [ dcoeffval		dDval ];	% ...and its error...
   coeffnames	= { coeffnames{:}	Dname };	% ...and its name.
  end
 end

end	% fit_discrete
