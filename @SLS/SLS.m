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
  C										% protein concentration
  Unit_C									% units
  Angles									% angles for SLS
  Q										% corresponding Q vectors
  Unit_Q									% units (A^{-1})
  Q2										% squared scattering vectors
  Unit_Q2									% units (A^{-2})
  Unit_KcR									% units ( Da^{-1})
  T										% temperature ( K )

  % N.B.: some dynamic properties could be called later on

 end	% public properties

 %============================================================================
 % PUBLIC PROPERTIES (READ-ONLY)
 %============================================================================
 properties ( SetAccess = private )

  % data struct
  Data

 end	% public properties

 %============================================================================
 % PUBLIC METHODS
 %============================================================================
 methods ( Access = public )

  %============================================================================
  % CONSTRUCTOR
  %============================================================================
  function sls = SLS ( dataclasses )

   % eat the common strings from the first Data class
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
   sls.Data.Angles	= horzcat(dataclasses.Angles_static);
   sls.Data.Q		= horzcat(dataclasses.Q_static);
   sls.Data.Q2		= horzcat(dataclasses.Q_static) .^2;
   sls.Data.KcR		= horzcat(dataclasses.KcR);
   % check for 0.00 errors in dKcR ( the ALV machine does it sometimes )
   min_relerr		= 1e-4; % mol/g	
   sls.Data.dKcR	= max(horzcat(dataclasses.dKcR),min_relerr*horzcat(dataclasses.KcR));

   % create unique properties with the update_uniques function in the Utils class
   sls.update_uniques('C','Q','Q2','Angles');

   % monitor properties with listeners in the Utils class
   sls.monitorprop('C','Q','Q2','Angles');

  end	% constructor

  %============================================================================
  % PLOT KC/R VERSUS CONCENTRATION, ANGLE, OR IONIC STRENGTH
  %============================================================================
  function plot_KcR ( self, varargin )
  % plot Kc/R in a parametric fashion
  %
  % if no optional argument is given, the x axis is the concentration, and
  % all angles are plotted together in one figure
  %
  % allowed arguments are the following:
  % - 'Independent',<x> where <x> can be 'C','I','Q2'
  % - 'Parameter', <p> where <p> can be 'C','I','Q2' (plot of C vs I not allowed)
  % - 'Color', <c> where <c> can be 'blue','green','red','pink'
  % - 'Figure',<f> where <f> must be a figure handle (e.g. gcf)

   options = struct(varargin{:});								% get the optional arguments

   try	independent	= options.Independent;	catch	independent	= 'C';	end		% default independent is C
   try	average		= options.Average;	catch	average	= 'yes';	end		% average if not prohibited
   try	color		= options.Color;	catch	color	= 'Green';	end		% default color is green

   try	fig		= options.Figure;							% try using existing figure...
   catch
    fig	= self.create_figure;									% ...or create one
    xlabel(self.set_xlabel(independent));							% choose x label
    ylabel(['Kc/R [ ',self.Unit_KcR,' ]']);							% y axis is always the scattering ratio
   end

   colos = struct('Green',[],'Red',[],'Blue',[],'Pink',[]);					% colors

    switch average

     case 'no'											% if no averaging is chosen...

      try	parameter	= options.Parameter;						% ...try to read parameter or...
      catch	parameter	= self.set_parameter(independent);				% ...fallback parameter
      end

      nuances	= self.get_color(color,length(self.(parameter)));				% get the colors for the plot

      for j = 1 : length(self.(parameter))
       inde	= self.Data.(parameter) == self.(parameter)(j);					% find the right data to plot
       x	= self.Data.(independent)	(inde);
       KcR	= self.Data.KcR			(inde);
       dKcR	= self.Data.dKcR		(inde);

       plo(j)	= errorbar( x, KcR, dKcR,	'o-','Color',nuances(j,:),		...
						'MarkerSize',	self.MarkerSize,	...
						'LineWidth',	self.LineWidth		);	% plot
       legend_names{j}	= [parameter,' = ',num2str(self.(parameter)(j),1),' ',...
								self.(['Unit_',parameter])];	% legend names
      end
      legend(plo,legend_names,'Location','NorthWest');						% plot legend

     case 'yes'											% if averaging is chosen...

      nuances	= self.get_color(color,1);							% get the colors for the plot

      for i = 1 : length(self.(independent))
       inde	= self.Data.(independent) == self.(independent)(i);				% find the right data to plot
       x(i)	= self.(independent)(i);
       KcR(i)	= mean ( self.Data.KcR ( inde ) );
       dKcR(i)	= mean ( self.Data.dKcR ( inde ) );						% INCORRECT!
      end

      plo	= errorbar( x, KcR, dKcR,	'o-','Color',nuances,			...
						'MarkerSize',	self.MarkerSize,	...
						'LineWidth',	self.LineWidth		);	% plot

    end

  end % plot_KcR
    
  %============================================================================
  % PLOT OSMOTIC COMPRESSIBILITY VERSUS CONCENTRATION OR IONIC STRENGTH
  %============================================================================
  function plot_compressibility ( self, varargin )
  % the osmotic isothermal compressibility is related to the structure factor at
  % vanishing Q by the following formula:
  %
  %	X := ( dc / dP )_T = N_A / [ M kT * S(c, q -> 0) ] = N_A / [ kT  (Kc/R) ]
  %
  % where X is the compressibility, c is the concentration in mass/volume, N_A is
  % Avogadro's number and S(c,q) is the static structure factor.
  %
  % Please note that X(c) is independent on Q and on M.

   options = struct(varargin{:});								% try to interpret input options

   try	independent	= options.Independent;	catch	independent	= 'C';	end		% default independent is C
   try	color		= options.Color;	catch	color	= 'Green';	end		% default color is green

   for i = 1 : length(self.(independent))
    KcR(i)	= mean ( self.Data.KcR ( self.Data.(independent) == self.(independent)(i) ) );	% mean of the scattering ratios at that concentration
    dKcR(i)	= mean ( self.Data.dKcR ( self.Data.(independent) == self.(independent)(i) ) );	% INCORRECT!
    X(i)	= LIT.Na / ( LIT.kb * self.T * KcR(i) );					% set compressibility
    dX(i)	= X(i) * dKcR(i) / KcR(i);							% set error
   end

   try fig	= options.Figure;								% try using existing figure...
   catch
    fig	= self.create_figure;									% ...or create one
    xlabel(self.set_xlabel(independent));							% choose the x label
    ylabel(['Compressibility [ Da J^{-1} ]']);							% choose the y label
   end

   nuances	= self.get_color(color,1);							% get the colors for the plot

   errorbar(self.(independent),X,dX,	'o-', 				...
					'Color',	nuances,	...
					'LineWidth',	self.LineWidth,	...
					'MarkerSize',	self.MarkerSize	);			% plot!

  end	% plot_compressibility

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
   coeff_units	= {'Da^{-1}'	'mol*ml/g^2'};
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
  % check whether a certain property is already in the class or not,
  % and add it if requested
  %============================================================================
  function check_add_prop ( obj, varargin )
   if mod(length(varargin),2)
    error('You should choose parameter names and input a values for them.');
   else
 
    props = struct(varargin{:});		% create the struct
    f = fieldnames(props);			% find out the fieldnames

    for i = 1 : length(f)
     if ~ismember(f{i},properties(obj))		% add the prop if necessary
      obj.addprop(f{i});
     end
     obj.(f{i}) = props.(f{i});			% fill the prop
    end

   end 
  end	% check_add_prop

 end	% public methods

 %============================================================================
 % PRIVATE METHODS
 %============================================================================
 methods ( Access = private )

  %============================================================================
  % set xlabel according to independent
  %============================================================================
  function xl = set_xlabel ( self, independent )

   switch independent
    case 'C'							% concentration
     xl = ['Protein concentration [ ',self.Unit_C,' ]'];

    case 'I'							% ionic strength
     xl = ['ionic strength [ ',self.Unit_I,' ]'];

    case 'Q2'							% Q squared
     xl = ['Q^2 [ ',self.Unit_Q2,' ]'];

    otherwise							% not recognized?!
     error('Independent not recognized!');
   end

  end	% set_xlabel

  %============================================================================
  % set parameter from independent
  %============================================================================
  function parameter = set_parameter ( self, independent )

   switch independent
    case 'C'							% concentration
     parameter = 'Q2';

    case 'I'							% ionic strength
     parameter = 'Q2';

    case 'Q2'							% Q squared
     parameter = 'C';

    otherwise							% not recognized?!
     error('Independent not recognized!');
   end

  end	% set_parameter

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
