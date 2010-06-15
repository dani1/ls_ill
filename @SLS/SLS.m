classdef SLS < dynamicprops & Graphics & Utils
%============================================================================
%
% SLS class for plotting, fitting, showing static light scattering data
% taken from one Data class or from a vector of classes.
%
% Author:	Fabio Zanini
% Companies:	Institut Laue Langevin, Universität Tübingen
% Date:		12 May 2010
% Version:	1.2
% Copyright:	2009,2010 Fabio Zanini
% License:	GPLv3
%
% Syntax: slsclass = SLS ( input_data_classes )
%
%
% A list of things you can do:
% - show data parameters (angles, q vectors, concentrations, etc.);
% - plot data in different forms (like at some angles or concentrations only);
% - fit the data according to usual SLS models.
%
%
% That is how Kc/R data are generated from the instrument. One measures
% essentially count rates. For these to be used, one needs to input following
% parameters into the SLS program:
%
% - sample concentration [mg/ml]
% - sample refractive index
% - sample refractive index increment with c [ml/g]
% - solvent type + viscosity
% - temperature
%
% Some of these parameters are unknown when one performs the measurement, so
% it is necessary to correct them later. This is *already* done in the Data
% class for every measurement individually, so no further change is needed.
% Anyway, the formula for calculating Kc/R would be:
%
%		Kc/R = ( 2 pi n dn/dc )^2 / ( lambda^4 * N_a )
%
% where lambda is the laser *vacuum* wavelength and Na is the Avogadro
% number.
%
%============================================================================

 %============================================================================
 % PUBLIC PROPERTIES
 %============================================================================
 properties ( Access = public )

  % general parameters
  Sample
  Unit_C									% units
  Unit_Q									% units (A^{-1})
  Unit_Q2									% units (A^{-2})
  Unit_KcR									% units ( Da^{-1})
  T										% temperature ( K )
  X_T										% osmotic compressibility
  dX_T										% error on X_T
  Unit_X_T									% units for compressibility
  Data										% data stored as fields (vectors or cell arrays)


  % N.B.: some dynamic properties could be called later on

 end	% public properties

 %============================================================================
 % OBSERVED PUBLIC PROPERTIES
 %============================================================================
 properties ( Access = public, SetObservable, AbortSet )
 % These are useful for triggering events after change in a property value
 % In this class, the Data struct is changed according to changes in Q and C

  C
  Phi
  Angles
  Q
  Q2

 end	% observed public properties

 %============================================================================
 % PUBLIC METHODS
 %============================================================================
 methods ( Access = public )

  %============================================================================
  % CONSTRUCTOR
  %============================================================================
  function sls = SLS ( dataclasses )

   % eat the common strings from the first Data class
   sls.Sample		= dataclasses(1).Sample;
   sls.Unit_KcR		= dataclasses(1).Unit_KcR;
   sls.Unit_C		= dataclasses(1).Unit_C;
   sls.Unit_Q		= dataclasses(1).Unit_Q;
   sls.Unit_Q2		= regexprep(sls.Unit_Q,'\^\{*-1\}*','\^\{-2\}');
   sls.T		= dataclasses(1).T;

   % here I try to introduce a structure similar to the class DLS, in order to
   % control the data flux through logical operators for index arrays, rather
   % than doing everything manually and turning the world upside down if one
   % wants to plot the data as a function of angle
   conc = [];
   for i = 1 : length(dataclasses)
    for j = 1 : length(dataclasses(i).Angles_static)
     conc	= [ conc dataclasses(i).C ];
    end
   end
   sls.Data		= struct('C',[],'Angles',[],'Q',[],'Q2',[],'KcR',[],'dKcR',[]);
   sls.Data.C		= conc;
   sls.Data.Phi		= sls.Data.C * LIT.(sls.Sample).v0;
   sls.Data.Angles	= horzcat(dataclasses.Angles_static);
   sls.Data.Q		= horzcat(dataclasses.Q_static);
   sls.Data.Q2		= horzcat(dataclasses.Q_static) .^2;
   sls.Data.KcR		= horzcat(dataclasses.KcR);
   % check for 0.00 errors in dKcR ( the ALV machine does it sometimes )
   min_relerr		= 1e-4; % mol/g	
   sls.Data.dKcR	= max(horzcat(dataclasses.dKcR),min_relerr*horzcat(dataclasses.KcR));

   [ sls.X_T sls.dX_T sls.Unit_X_T ] = sls.calculate_compressibility;		% calculate the compressibility

   % create unique properties with the update_uniques function in the Utils class
   sls.update_uniques('C','Phi','Q','Q2','Angles');

   % monitor properties with listeners in the Utils class
   sls.monitorprop('C','Phi','Q','Q2','Angles');

  end	% constructor

  %============================================================================
  % CALCULATE OSMOTIC COMPRESSIBILITY
  %============================================================================
  function [ X dX Unit_X ] = calculate_compressibility ( self )
  % the osmotic isothermal compressibility is related to the structure factor at
  % vanishing Q by the following formula:
  %
  %	X :=	1/c * ( dc / dP )_T	= 1 / [ Na * kT * c * M * S(c, q -> 0)	]
  %					= 1 / [ Na * kT * c * 	(Kc/R) 		]
  %
  % where X is the compressibility, c is the concentration in mass/volume, N_A is
  % Avogadro's number and S(c,q) is the static structure factor.
  %
  % Please note that X(c) is independent on Q and on M.

   Na	= LIT.Constants.Na;
   kb	= LIT.Constants.kb;
   T	= self.T;
   C	= self.Data.C;
   KcR	= self.Data.KcR;
   dKcR	= self.Data.dKcR;

   X	= 1 ./ ( Na .* kb .* T .* C .* KcR );			% calculate X_T
   dX	= X .* dKcR ./ KcR;					% propagate errors
   Unit_X = 'l * J^{-1}';

  end	% calculate_compressibility

  %============================================================================
  % FIT SCATTERING RATIO VS CONCENTRATION (A2) ETC.
  %============================================================================
  function fit_KcR ( obj , varargin )
  % fit static data through M, Rg
  % Input values:
  % - class (SLS)
  % - optional arguments for fine-tuning.
  %
  % Optional arguments are the following:
  % - 'N', n:			the first n points are taken for the fit (default: all)
  % - 'Plot', 'yes'/'no':	plot the fit afterwards
  %
  % One gets the fit parameters M, B2 and their fit uncertainties.
  %
  % Explanation: the mean scattered intensity should follow the law
  %
  %		Kc/R = 1/M * ( 1 + 0.33 * Rg^2 q^2 ) + 2 * B2 * c
  %
  % Therefore, one should find the mass from the zero-q and zero-c
  % extrapolation, and the other two parameters from the complete fit.
  %
  % Obviously, the small-q and small-c limits *should* be order-independent...

   % understand which parameters are selected
   p = struct(varargin{:});
   p = obj.check_fit( p );   

   % decide how many points are taken for the linear fit
   if ~isfield(p,'N')
      p = setfield(p,'N',length(obj.C));
   end

   datax	= obj.C(1:p.N);							% create the x dataset from the number of points

   for i = 1 : length(datax)
    inde	= obj.Data.C == datax(i);					% find the right data to fit
    datay(i)	= mean ( obj.Data.KcR ( inde ) );
    datady(i)	= mean ( obj.Data.dKcR ( inde ) );
   end

   fit_function	= 'Minv + 2 * B2 * C';						% fit function
   coefficients	= {'Minv'	'B2'};						% coefficients
   coeff_units	= {'Da^{-1}'	'mol*l/g^2'};
   fit_expr	= {'1'		'2*C'};						% one MUST use this syntax with linear models

   weights	= 1 ./ ( datady ) .^2;

   fit_options	= fitoptions('Method','LinearLeastSquares','Weights', weights);	% fit options
   fit_type	= fittype(fit_expr, 'independent','C',...
		'coefficients',coefficients,'options',fit_options);		% fit type

   [ cf gof ] = fit( datax', datay', fit_type );				% perform numerical fit

   % get fit result
   coeffout = coeffvalues (cf); 
   confidence = confint(cf,0.95);
   dcoeffout =  0.5*(confidence(2,:)-confidence(1,:));  

   % output coefficients: mass
   M	= 1 / coeffout(1);							% mass
   dM	= dcoeffout(1) / coeffout(1) * M;					% mass error
   obj.check_add_prop('M',M,'dM',dM,'Unit_M','Da');				% mass unit

   for k = 2 : length(coefficients)
    obj.check_add_prop(	coefficients{k},coeffout(k),			...
			['d',coefficients{k}],dcoeffout(k),		...
			['Unit_',coefficients{k}],coeff_units{k}	);	% output other coefficients (A2, etc.) to the class
   end

    fit_annot	= { ['M = ',num2str(obj.M),' ',obj.Unit_M],		...
		[coefficients{2},' = ',num2str(coeffout(2)),		...
		 ' ', coeff_units{2}] };					% fit_annotation (for plot)

   % plot the fit if requested
   if strcmpi(p.Plot,'yes')

    obj.plot_KcR('Independent','C','Average','yes');				% plot the data

    fun_plot	= inline( fit_function, 'C', coefficients{:}	);		% create the inline function for the fit

    colos	= get(get(gca,'Children'),'Color');
    colos	= colos + 0.5*(1-colos);

    plot( obj.C, 						...
	fun_plot( obj.C,coeffout(1),coeffout(2) )	,	...
	'--', 'Color',colos,'LineWidth',2				);	% plot the fit

    text( min(datax), max(datay), fit_annot, 'FontSize',14);			% text annotation
   end
   
  end	% fit_KcR

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

  %==========================================================================
  % ERRORBAR (OVERLOADED METHOD)
  %==========================================================================
  function errorbar ( self, name, varargin )
  % function for plotting every kind of property quickly
  % first we look in the Data struct. Otherwise, the property is searched for directly
  % Accepted optional arguments:
  % - Independent:	name of the independent
  % - Average:		averaging over some kind of parameter?
  % - Parameter:	what is the parameter ( in case of parametric plots)
  % - Color:		what color nuances to use
  %
  % N.B.: this function can be overloaded by subclass methods, and should be considered
  % a fallback method

   optargs	= varargin;								% get the optional arguments

   [ ax, x, y, dy, color ] = self.prepare_errorbar ( name, 'Type', 'errorbar', optargs{:} );	% prepare the plot

   if strcmp(name,'X_T')     set(ax,'XScale','log','YScale','log');		end	% set loglog for compressibility

   self.make_errorbar ( ax, x, y, dy, color, optargs{:} );					% make the plot

  end % errorbar


  %==========================================================================
  % PLOT (OVERLOADED METHOD)
  %==========================================================================
  function plot ( self, name, varargin )
  % function for plotting every kind of property quickly
  % first we look in the Data struct. Otherwise, the property is searched for directly
  % Accepted optional arguments:
  % - Independent:	name of the independent
  % - Average:		averaging over some kind of parameter?
  % - Parameter:	what is the parameter ( in case of parametric plots)
  % - Color:		what color nuances to use
  %
  % N.B.: this function can be overloaded by subclass methods, and should be considered
  % a fallback method

   optargs	= varargin;								% get the optional arguments

   switch name

    case 'X_T'
     [ ax, x, y, color ] = self.prepare_plot ( name, optargs{:} );			% prepare the plot
     set(ax,'XScale','log','YScale','log');						% set loglog for compressibility
     self.make_plot ( ax, x, y, color, optargs{:} );						% make the plot

    case 'Fit_KcR'
     [ ax, x, y, color ] = self.prepare_plot ( 'KcR', optargs{:} );			% prepare the plot
     self.make_plot ( ax, x, y, color, optargs{:} );						% make the plot

     options	= struct(optargs{:});
     try independent = options.Independent;	catch independent = 'C'; end
     q	= 1 / self.M;									% select the right intercept
     m	= 2 * self.B2;									% select the right slope

     if strcmp(independent,'Phi') m = m/LIT.(Sample).v0;	end			% change the slope if necessary

     self.plot_affine_fit ( x, q, m );							% plot the fit

    

    otherwise
     [ ax, x, y, color ] = self.prepare_plot ( name, optargs{:} );			% prepare the plot
     self.make_plot ( ax, x, y, color, optargs{:} );					% make the plot

   end

  end % plot


 end	% public methods

 %============================================================================
 % PRIVATE METHODS
 %============================================================================
 methods ( Access = private )

  %============================================================================
  % check the parameters for fits
  %============================================================================
  function par = check_fit( obj, par )

   % choose the independent, the x axis (default: C)
   if ~isfield(par,'Independent')
    par = setfield(par,'Independent','C');
   end

   if strcmp('Independent','C')
    par = setfield(par,'Parameter','Q2')
   end

   % choose to plot or not the fit results
   if ~isfield(par,'Plot')
    par = setfield(par,'Plot','no');
   end

  end	% check_fit

 end	% private methods

end	% SLS
