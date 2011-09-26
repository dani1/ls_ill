classdef Experiment < dynamicprops
%============================================================================
%
% Experiment class for plotting, fitting, showing static light scattering data
% taken from one Data class or from a vector of classes.
%
% Author:	Fabio Zanini
% Companies:	Institut Laue Langevin, Universität Tübingen
% Date:		12 May 2010
% Version:	2.0
% Copyright:	2009,2010 Fabio Zanini
% License:	GPLv3
%
% Syntax: slsclass = Experiment ( vector_of_points )
%============================================================================

 %============================================================================
 % PUBLIC PROPERTIES
 %============================================================================
 properties ( Access = public )

  % general parameters
  Protein
  Unit_C									% units
  Unit_Cs									% units ( mM )
  Unit_Q									% units (A^{-1})
  Unit_Q2									% units (A^{-2})
  Unit_KcR									% units ( Da^{-1})
  T										% temperature ( K )
  Unit_T									% units for compressibility
  Unit_X_T									% units for compressibility
  Point										% data stored as fields (vectors or cell arrays)


  % N.B.: some dynamic properties could be called later on

 end	% public properties

 %============================================================================
 % DEPENDENT PROPERTIES
 %============================================================================
 properties ( GetAccess = public, Dependent, SetAccess = private )
 % These are useful for triggering events after change in a property value
 % In this class, the Data struct is changed according to changes in Q and C

  C
  Phi
  Cs
  Angle
  Q
  Q2

 end

 %============================================================================
 % PUBLIC METHODS
 %============================================================================
 methods ( Access = public )

  %============================================================================
  % CONSTRUCTOR
  %============================================================================
  function self	= Experiment ( classin )

   if strfind(class(classin),'Sample')						% if you call the class from samples, add both Point and Sample
    point	= [classin.Point];
    self.addprop('Sample');
    self.Sample	= [classin];
    self.addprop('KcR');
    self.KcR	= [ self.Sample.KcR ];
    self.addprop('dKcR');
    self.dKcR	= [ self.Sample.dKcR ];
   else										% otherwise, only add Point
    point	= [classin];
   end

   % eat the common strings from the first Data class
   self.Protein		= point(1).Protein;
   self.Unit_KcR	= point(1).Unit_KcR;
   self.Unit_C		= point(1).Unit_C;
   self.Unit_Cs		= point(1).Unit_Cs;
   self.Unit_Q		= point(1).Unit_Q;
   self.Unit_Q2		= regexprep(self.Unit_Q,'\^\{*-1\}*','\^\{-2\}');
   self.Unit_T		= point(1).Unit_T;
   self.Unit_X_T	= point(1).Unit_X_T;

   self.Point		= point;						% import the point classes as props

  end	% constructor

 end

 methods

  %============================================================================
  % GET METHODS FOR DEPENDENT PROPERTIES
  %============================================================================
  function C = get.C ( self )
   C	= unique( [ self.Point.C ] );
  end

  function Phi = get.Phi ( self )
   Phi	= unique( [ self.Point.Phi ] );
  end

  function Cs = get.Cs ( self )
   Cs	= unique( [ self.Point.Cs ] );
  end

  function Angle = get.Angle ( self )
   Angle	= unique( [ self.Point.Angle ] );
  end

  function Q = get.Q ( self )
   Q	= unique( [ self.Point.Q ] );
  end

  function T = get.T ( self )
   T	= mean( [ self.Point.T ] );
  end

 end

 methods
  %============================================================================
  % FIT SCATTERING RATIO VS CONCENTRATION (A2) ETC.
  %============================================================================
  function [ M B2 dM dB2 ] = fit_KcR ( self , N )
  % fit static data linearly, to find M and B2
  %
  % Optional arguments are the following:
  % - N:			the first n points are taken for the fit (default: min(4,all))
  %
  % Explanation: the mean scattered intensity should follow the law
  %
  %		Kc/R = 1/M + 2 * B2 * c
  %

   if nargin < 2
    N	= min(4,length(self.C));
   end

   index	= ( [ self.Point.C ] <= self.C(N) );
   point	= self.Point( index );

   C	= [ point.C	]';
   KcR	= [ point.KcR	]';
   dKcR	= [ point.dKcR	]';

   fitfun	= '1 ./ M + 2 .* B2 .* C';					% fit function
   coefficients	= {'Mi'	'B2'};						% coefficients
   fitexpr	= {'1'	'2*C'};

   weights	= 1 ./ ( dKcR ) .^2;

   fit_options	= fitoptions(	'Method',	'LinearLeastSquares',	...
				'Weights',	weights			);	% fit options

   fit_type	= fittype(	fitexpr, 				...
				'independent',	'C',			...
				'dependent',	'KcR',			...
				'coefficients',	coefficients,		...
				'options',	fit_options		);	% fit type

   [ cf gof ] = fit( C, KcR, fit_type );					% perform numerical fit
   confidence = confint(cf,0.95);

   % get fit result
   cout = coeffvalues (cf); 
   dcout =  0.5*(confidence(2,:)-confidence(1,:));  
   M	= 1/ cout(1);							% mass
   dM	= dcout(1) / cout(1)^2 ;
   B2	= cout(2);							% second virial coefficient
   dB2	= dcout(2);

   try	self.addprop('Fit_KcR');	end
   self.Fit_KcR	= cf;

  end

  %============================================================================
  % FIT COMPRESSIBILITY
  %============================================================================
  function fit_compressibility ( self, varargin )
  % fit the compressibility according to some nonlinear model

   options = struct(varargin{:});						% eat the optional arguments as a struct

   try N = options.N;	catch N = length(self.C); end				% how many points to use for the fit

   try
    x	= self.C	(1:N);								% x vector
    y	= self.X_T	(1:N);								% y vector
    dy	= self.dX_T	(1:N);								% dy vector
   catch
    error('Does X_T exist?');
   end

   weights =  1 ./ dy' .^ 2;							% set the weights for the fit

   fit_function	= 'BKG + beta * C ^ ( - eta )';
   coefficients	= {'BKG'	'beta'	'eta'	};
   lowerbond		= [ 0		0	0	];			% eta is between 0 (no C dependence)...
   upperbond		= [ +inf	+inf	2	];			% ...and 2 (finite potential at r-->oo)
   startpoint		= [ 1		1	1.5	];

   % set other fit options
   fit_options = fitoptions( 'Method','NonLinearLeastSquares',...
      		'Weights', weights,			...
      		'MaxFunEvals',100000, 			...
      		'Lower',lowerbond, 			...
      		'Upper',upperbond, 			...
      		'Startpoint',startpoint				);
   
   fit_type = fittype(fit_function, 'independent','C', ...
      		'options', fit_options,			...
      		'coefficients',	coefficients			);

   % perform the numerical fit
   [ cf gof ] = fit( x', y', fit_type );

   % get fit result
   par = coeffvalues (cf); 
   confidence = confint(cf,0.95);
   dpar =  0.5*(confidence(2,:)-confidence(1,:));  

   % export the fit values
   BKG		= par(1);
   beta	= par(2);
   eta		= par(3)
   deta	= dpar(3)

   self.plot_compressibility;
   plot(x, BKG + beta .* x.^( -eta ), 'LineWidth', self.LineWidth );

   self.check_add_prop('Eta',eta,'dEta',deta);

  end	% fit_compressibility

 end

end
