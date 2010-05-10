classdef SLS < dynamicprops
%============================================================================
%
% SLS class for plotting, fitting, showing static light scattering data
% taken from one LSData class or from a *cell array* of classes.
%
% Author:	Fabio Zanini
% Companies:	Institut Laue Langevin, Universität Tübingen
% Date:		22 Mar 2010
% Version:	1.1
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

   % eat concentrations, q, angles
   sls.C		= unique(sls.Data.C);
   sls.Angles		= unique(sls.Data.Angles);
   sls.Q		= unique(sls.Data.Q);
   sls.Q2		= sls.Q .^2;

  end	% constructor

  %============================================================================
  % PLOT KC/R VERSUS CONCENTRATION, ANGLE, OR IONIC STRENGTH
  %============================================================================
  function plot_KcR ( obj, varargin )
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

   % understand which parameters are selected, create the figure if requested
   p = struct(varargin{:});
   p = obj.check_plot( p );

   % create colors
   colos = struct('Green',[],'Red',[],'Blue',[],'Pink',[]);
   for j = 1 : length(obj.(p.Parameter))
    colos.Green(j,:)	= [ 0.1 j/length(obj.(p.Parameter)) 0.2 ];
    colos.Red(j,:)	= [ 0.8 j/length(obj.(p.Parameter)) 0.2 ];
    colos.Blue(j,:)	= [ 0.1 j/length(obj.(p.Parameter)) 0.8 ];
    colos.Pink(j,:)	= [ 0.7 j/length(obj.(p.Parameter)) 0.6 ];
   end

    switch p.Average
     case 'no'
      for j = 1 : length(obj.(p.Parameter))
       colo		= colos.(p.Color);									% color stuff
       inde		= obj.Data.(p.Parameter) == obj.(p.Parameter)(j);					% find the right data to plot
       plo(j)		= errorbar(obj.Data.(p.Independent)(inde),obj.Data.KcR(inde),obj.Data.dKcR(inde),...
							'o-','Color',colo(j,:),'MarkerSize',10,'LineWidth',4);	% plot
       legend_names{j}	= [p.Parameter,' = ',num2str(obj.(p.Parameter)(j),1),' ',...
										obj.(['Unit_',p.Parameter])];	% legend names
      end
      legend(plo,legend_names,'Location','NorthWest');								% plot legend

     case 'yes'
      for i = 1 : length(obj.(p.Independent))
       inde		= obj.Data.(p.Independent) == obj.(p.Independent)(i);					% find the right data to plot
       x(i)		= obj.(p.Independent)(i);
       KcR(i)		= mean ( obj.Data.KcR ( inde ) );
       dKcR(i)		= mean ( obj.Data.dKcR ( inde ) );
      end
      colo		= colos.(p.Color)(1,:);									% color stuff
      plo		= errorbar(x,KcR,dKcR,'o-','Color',colo,'MarkerSize',10,'LineWidth',4);			% plot

    end

  end % plot_KcR
    
  %============================================================================
  % PLOT OSMOTIC COMPRESSIBILITY VERSUS CONCENTRATION OR IONIC STRENGTH
  %============================================================================
  function plot_compressibility ( obj, varargin )
  % the osmotic isothermal compressibility is related to the structure factor at
  % vanishing Q by the following formula:
  %
  %	X := ( dc / dP )_T = N_A / [ M kT * S(c, q -> 0) ] = N_A / [ kT  (Kc/R) ]
  %
  % where X is the compressibility, c is the concentration in mass/volume, N_A is
  % Avogadro's number and S(c,q) is the static structure factor.
  %
  % Please note that X(c) is independent on Q and on M.

   try options = struct(varargin{:}); end									% try to interpret input options

   try	independent	= options.Independent;	catch	independent	= 'C';	end				% default independent is C

   for i = 1 : length(obj.(independent))
    KcR(i)	= mean ( obj.Data.KcR ( obj.Data.(independent) == obj.(independent)(i) ) );			% mean of the scattering ratios at that concentration
    dKcR(i)	= mean ( obj.Data.dKcR ( obj.Data.(independent) == obj.(independent)(i) ) );			% INCORRECT!
    X(i)	= LIT.Na / ( LIT.kb * obj.T * KcR(i) );								% set compressibility
    dX(i)	= X(i) * dKcR(i) / KcR(i);									% set error
   end

   try fig	= options.Figure;										% try using existing figure...
   catch
    fig	= obj.fig_KcR_c;											% ...or create one
    ylabel(['Compressibility [ Da J^{-1} ]']);								% modify the label if creating a new figure
   end

   errorbar(obj.(independent),X,dX,'o-','LineWidth',3,'MarkerSize',10);						% plot!

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
  % create a figure for Kc/R
  %============================================================================
  function fig = fig_KcR_c ( obj )
   
   fig = figure;
   hold all;
   set( gcf, 'color', 'white' );
   set( gca,	'box', 'on', 						...
		'FontSize',36,'FontName','Courier','FontWeight','bold',	...
		'XMinorTick','on', 'YMinorTick','on');
   ylabel(['Kc/R [ ',obj.Unit_KcR,' ]']);

  end	% fig_KcR_c
  
  %============================================================================
  % check the parameters for plots
  %============================================================================
  function par = check_plot( obj, par )

   % set the figure and add the property to the class
   if ~isfield(par,'Figure')
    fig = obj.fig_KcR_c;
    par=setfield(par,'Figure',fig);
    obj.check_add_prop('Figure',par.Figure);
   end

   % choose the independent, the x axis (default: C)
   if ~isfield(par,'Independent')
    par = setfield(par,'Independent','C');
   end
   switch par.Independent
    case 'C'
     xlabel(['Protein concentration [ ',obj.Unit_C,' ]']);
     if ~isfield(par,'Parameter')
      par = setfield(par,'Parameter','Q2');
     end
     if ~isfield(par,'Average')
      par = setfield(par,'Average','no');
     end
    case 'I'
     xlabel(['ionic strength [ ',obj.Unit_I,' ]']);
     par = setfield(par,'Parameter','Q2');
    case 'Q2'
     xlabel(['Q^2 [ ',obj.Unit_Q2,' ]']);
     if ~isfield(par,'Parameter')
      par = setfield(par,'Parameter','C');
     end
     if ~isfield(par,'Average')
      par = setfield(par,'Average','no');
     end
    otherwise
     error('Independent not recognized!');
   end

   % title
   title(['Scattering ratio as a function of ',par.Independent]);

   % here various colors
   if ~isfield(par,'Color') 
    par = setfield(par,'Color','Green');
   end

  end	% check_plot

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
