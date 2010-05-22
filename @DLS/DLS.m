classdef DLS < dynamicprops & Graphics & Utils
%============================================================================
%
% DLS class for plotting, fitting, showing dynamic light scattering data
% taken from an array of type Data
%
% Author:	Fabio Zanini
% Companies:	Institut Laue Langevin, Universität Tübingen
% Date:		15 May 2010
% Version:	1.3
%
% Syntax: dlsclass = DLS( input_data_classes )
%
% A list of things you can do:
% - show data parameters (angles, q vectors, concentrations, etc.);
% - plot data;
% - fit the data according to usual DLS models (single or double exponential,
%   cumulant of second order, streched exponential);
% - find the H function using both DLS and SLS.
%
%============================================================================

 %============================================================================
 % PUBLIC PROPERTIES
 %============================================================================
 properties (Access=public)
 % N.B.: some dynamic properties could be called later on, like Gamma, D, etc.

  Unit_C
  Unit_Q
  Unit_Q2
  % both Angles and Q are the unique reduction (union) of all their partners
  % in the input Data classes
  Unit_Tau
  T

 end	% public properties
 
 %============================================================================
 % OBSERVED PUBLIC PROPERTIES
 %============================================================================
 properties ( Access = public, SetObservable, AbortSet )
 % These are useful for triggering events after change in a property value
 % In this class, the Data struct is changed according to changes in Q and C

  C
  Angles
  Q
  Q2

 end	% observed public properties

 %============================================================================
 % HIDDEN PUBLIC PROPERTIES
 %============================================================================
 properties ( Access = public, Hidden)
  Data
  % Data is a struct with all the data stored as fields (vectors or cell arrays)
  % For instance, all correlations are stored in dls.Data.G
  % Therefore one is allowed to use such syntaxes like the following:
  %
  % dls.Data.G ( dls.Data.Angles == 30 )
  %
  % meaning that we look for all correlogramms for data at 30°. This is *very*
  % flexible and the best way to store things (remember dynamic field get-set
  % methods, etc.)
  %
  Minimal_D1		= 1		% A^2 / ns
  Maximal_D1		= 100		% A^2 / ns
  Minimal_D2		= 0		% A^2 / ns
  Maximal_D2		= 2.5		% A^2 / ns

  MinimalLagtime	= 1e-6		% ms
  MaximalLagtime	= 1e4		% ms
  MinimalG		= -0.05
  MaximalG		= 1.2

 end	% properties

 %============================================================================
 % PUBLIC METHODS
 %============================================================================
 methods ( Access = public )

  %============================================================================
  % CONSTRUCTOR
  %============================================================================
  function dls = DLS ( dataclasses )
  % this class takes an array of Data classes as input

   % eat the strings
   dls.Unit_Q	= dataclasses(1).Unit_Q;
   dls.Unit_Q2	= regexprep(dls.Unit_Q,'\^\{*-1\}*','\^\{-2\}');
   dls.Unit_C	= dataclasses(1).Unit_C;
   dls.Unit_Tau	= dataclasses(1).Unit_Tau;
   dls.T	= dataclasses(1).T;

   % write the data as a structure (easier to manage)
   conc = [];
   for i = 1 : length(dataclasses)
    for j = 1 : length(dataclasses(i).Angles_dyn)
     conc	= [ conc dataclasses(i).C ];
    end
   end
   dls.Data		= struct('C',[],'Angles',[],'Q',[],'Tau',[],'G',[],'dG',[]);
   dls.Data.C		= conc;
   dls.Data.Angles	= horzcat(dataclasses.Angles_dyn);
   dls.Data.Q		= horzcat(dataclasses.Q_dyn);
   dls.Data.Q2		= horzcat(dataclasses.Q_dyn).^2;
   dls.Data.Tau		= horzcat(dataclasses.Tau);
   dls.Data.G		= horzcat(dataclasses.G);
   dls.Data.dG		= horzcat(dataclasses.dG);

   % create unique properties with the update_uniques function in the Utils class
   dls.update_uniques('C','Q','Q2','Angles');

   % monitor properties with listeners in the Utils class
   dls.monitorprop('C','Q','Q2','Angles');

  end % constructor

  %============================================================================
  % PLOT CORRELATIONS
  %============================================================================
  function plot_correlations ( self, varargin )
  % plot the autocorrelation data
  %
  % varargin allows the user to choose some details of the plot. The syntax is
  % e.g. the following:
  %
  % plot_correlations('figure',gcf,'C',[ 2 3 ],'Angles',[30 50 100])
  %
  % The optional parameters are the following (default values in brackets):
  % - figure (new figures)
  % - Angles (all angles)
  % - C (all concentrations)

   options = struct(varargin{:});						% eat the optional arguments as a struct

   try	C = options.C;			catch	C = self.C;		end	% concentrations
   try	angles = options.Angles;	catch	angles = self.Angles;	end	% angles

   try	color = options.Color;		catch   color	= 'Green';	end	% color
   nuances	= self.get_color(color,length(C));				% nuances for plot colors

   
   % create the figures and plot the data
   for j = 1 : length(angles)							% for every angle...

    % create the figure
    fig(j)		= self.create_figure;					% ...create the figure...
    ax(j)		= get(fig(j),'CurrentAxes');				% ...get the axis...

    set(gca, 'xscale','log');							% ...set the log scale...
    axis([ 	self.MinimalLagtime	self.MaximalLagtime	...
		self.MinimalG		self.MaximalG 		]);		% ...set the limits for the plot...

    xlabel(['\tau [ ',self.Unit_Tau,' ]']);					% ...set the x label...
    ylabel('g_1 [ normalised ]');						% ...set the y label...
    title(['Correlogramms at ',num2str(angles(j)),'°']);			% ...set a title...

    legend_names{j}	= {};	legend_group{j}	= [];				% ...and prepare legends

    for l = 1 : length(self.Data.Tau)						% look for data to be plotted...

     if ( self.Data.Angles(l) == angles(j) && ismember(self.Data.C(l),C) )	% ...look the right indices
 
      x		= self.Data.Tau{l};						% ...rename the data
      y		= self.Data.G{l};
      dy	= self.Data.dG{l};

      nuance	= nuances( self.Data.C(l) == C, :);				% select the nuance for the plot

      plo	= errorbar(ax(j), x, y, dy, 		'-o',		...
					'Color',	nuance,		...
					'LineWidth', 0.1*self.LineWidth,...
					'MarkerSize',0.3*self.MarkerSize );	% plot

      legen	=['Protein conc = ',num2str(self.Data.C(l),1),' ',self.Unit_C];	% legend

      % create legend arrays
      if ~ismember({legen},legend_names{j})
       legend_names{j}	= [ legend_names{j}; {legen}	];
       legend_group{j}	= [ legend_group{j}; plo	];
      end

     end
    end

    legend(get(fig(j),'CurrentAxes'),legend_group{j},legend_names{j},'FontSize',14);	% show the legend for every figure

   end

  end	% end of plot function

  %============================================================================
  % FIT CORRELATIONS
  %============================================================================
  function fit_correlations ( obj , method, varargin )
  % This functions fits all the correlogramms of the class. Method is one of the following:
  % - 'Single': single exponential;
  % - 'Double'(default): double exponential;
  % - 'Single+Streched': single fast exponential and slow streched exponential (0<beta<1);
  % Not in use but easy to implement:
  % - 'Cumulant2', cumulant second order (see Koppel et al.);
  % - 'Cumulant3', cumulant third order.
  %
  % accepted optional arguments: 'fit_gamma_q2' ('yes' or 'no')

   if nargin == 1	method = 'Double';   end		% default method: double exponential
 
   options = struct(varargin{:});				% interpret the optional arguments

   for l = 1 : length(obj.Data.Tau)				% fit for all data
 
     % select which range of data to use for the fit. The problem is mainly
     % that for long lag times I'm fitting background, whereas for very short
     % lag times it is afterpulse noise of the APDs
     inde	= obj.Data.Tau{l} > 50e-6 & obj.Data.Tau{l} < 1e2;
     tau	= obj.Data.Tau{l}(inde);
     g		= obj.Data.G{l}(inde);
     dg		= obj.Data.dG{l}(inde);

     weights =  1 ./ dg .^ 2;					% set the weights for the fit

     [	fit_function,	...
        coefficients,	...
        lowerbond,	...
	upperbond,	...
	startpoint  	] = obj.fit_parameters( method, l );	% find the limits for the fit parameter

     % set other fit options
     fit_options = fitoptions( 'Method','NonLinearLeastSquares',...
			'Weights', weights,			...
			'MaxFunEvals',100000, 			...
        		'Lower',lowerbond, 			...
			'Upper',upperbond, 			...
			'Startpoint',startpoint				);
     
     fit_type = fittype(fit_function, 'independent','t', 	...
			'options', fit_options,			...
        		'coefficients',	coefficients			);

     % perform the numerical fit
     [ cf gof ] = fit( tau, g, fit_type );

     % get fit result
     par = coeffvalues (cf); 
     confidence = confint(cf,0.95);
     dpar =  0.5*(confidence(2,:)-confidence(1,:));  

     % output all coefficients using dynamic variable names
     for k = 1 : length(coefficients) 
      coout(k,l)	= par(k);
      dcoout(k,l)	= dpar(k);
     end

   end

   % TODO: control the behaviour of the correlations at 100 times higher than the
   %        first decay time. If there is still something important, do not fit
   %        the data and impose NaN for both Gamma and dGamma. One can not live with
   %        scorpions is his bed!


   obj.check_add_prop('Unit_Gamma',	[ obj.Unit_Tau,'^{-1}']);			% add Unit_Tau to the class
   obj.check_add_prop('Unit_D',		'A^2/ns');					% add Unit_D to the class
   for k = 1 : length(coefficients)							% add other coefficients
    obj.check_add_prop(coefficients{k},coout(k,:));
    obj.check_add_prop(['d',coefficients{k}],dcoout(k,:));

    if regexp(coefficients{k},'Gamma')							% add diffusion constants...
     obj.check_add_prop(['D',coefficients{k}(6),'_a'], ...
					1e-6 .* coout(k,:) ./ obj.Data.Q.^2 );		% ...as Gamma/Q^2 and...
     obj.check_add_prop(['dD',coefficients{k}(6),'_a'], ...
					1e-6 .* dcoout(k,:) ./ obj.Data.Q.^2 );		% ...as dGamma/Q^2
    end
   end

   try	fitgammas = options.Fit_Gamma_Q2;   catch	fitgammas = 'no';   end		% default: do not perform the gamma fit
   if strcmp(fitgammas,'yes');								% if requested, fit the gammas

     % fit all gammas
     for k = 1 : length(coefficients)
      if regexp(coefficients{k},'Gamma')
       obj.fit_gamma_q2(coefficients{k});
      end
     end

   end

  end	% fit_correlations

  %============================================================================
  % FIND LIMITS FOR CORRELATION FIT
  %============================================================================
  function [ fit_function coefficients lowerbond upperbond startpoint ] = fit_parameters ( self, method, i )
  % set the lower and upper limits for the correletion fit according to Q^2 and the class props
  % the startpoint is taken as D0 for the faster D, and 0.1*D0 for the slower one

   % set the fit parameters, according to the fit method chosen
   % general parameters (background)
   fit_common	= 'BKG';
   coeff_common	= {'BKG'};
   lower_common	= 0;
   upper_common	= 1e-1;
   start_common	= 1e-3;

   switch method
 
    case 'Single'
     fit_function	= [ fit_common,'+ Ae * exp( - 2 * Gammae * t )'			];
     coefficients	= [ coeff_common	{'Ae'}	{'Gammae'}			];
     lowerbond		= [ lower_common	0.8	1e6*self.Minimal_D1*self.Data.Q2(i)	];
     upperbond		= [ upper_common	1.2	1e6*self.Maximal_D1*self.Data.Q2(i)	];
     startpoint		= [ start_common	1.0	1e6*LIT.BSA.D0*self.Data.Q2(i)		];

     case 'Streched'
     fit_function	= [ fit_common,'+ As * exp( - 2 * (Gammas * t)^b )'	];
     coefficients	= [ coeff_common	{'As'}	{'Gammas'}				{'b'}	];
     lowerbond		= [ lower_common	0.8	1e6*self.Minimal_D1*self.Data.Q2(i)	0	];
     upperbond		= [ upper_common	1.2	1e6*self.Maximal_D1*self.Data.Q2(i)	1	];
     startpoint		= [ start_common	1.1	1e6*LIT.BSA.D0*self.Data.Q2(i)		0.8	];
 
    case 'Double'
     fit_function	= [ fit_common,'+ ( A1 * exp( - Gamma1 * t ) + A2 * exp( - Gamma2 * t ) ).^2'	];
     coefficients	= [ coeff_common	{'A1'}	{'Gamma1'}				{'A2'}	{'Gamma2'}				];
     lowerbond		= [ lower_common	0.5	1e6*self.Minimal_D1*self.Data.Q2(i)	0	1e6*self.Minimal_D2*self.Data.Q2(i)	];
     upperbond		= [ upper_common	1	1e6*self.Maximal_D1*self.Data.Q2(i)	0.5	1e6*self.Maximal_D2*self.Data.Q2(i) 	];
     startpoint		= [ start_common	0.9	1e6*LIT.BSA.D0*self.Data.Q2(i)		0.1	1e5*LIT.BSA.D0*self.Data.Q2(i)		];
 
    case 'Single+Streched'
     fit_function	= [ fit_common,'+ ( Ae * exp( - Gammae * t ) + As * exp(-( Gammas *t).^b )) .^2'];
     coefficients	= [ coeff_common	{'Ae'}	{'Gammae'}				{'As'}	{'Gammas'}				{'b'}	];
     lowerbond		= [ lower_common	0.5	1e6*self.Minimal_D1*self.Data.Q2(i)	0	1e6*self.Minimal_D2*self.Data.Q2(i)	0	];
     upperbond		= [ upper_common	1	1e6*self.Maximal_D1*self.Data.Q2(i)	0.5	1e6*self.Maximal_D2*self.Data.Q2(i)	1	];
     startpoint		= [ start_common	0.9	1e6*LIT.BSA.D0*self.Data.Q2(i)		0.1	1e5*LIT.BSA.D0*self.Data.Q2(i)		0.8	];

    otherwise
     error('Method not recognized!');
   end

  end	% find_limits_for_fit
  
  %============================================================================
  % FIT CUMULANT (DEAD)
  %============================================================================
  % function varargout = fit_cumulant ( obj, order )
  % fit autocorrelation data according to Koppel JCP 57,11, 4814 (1972)
  %
  % the fit tries to fit with a cumulant expansion, plus amplitude
  % normalization, plus background, i.e. the assumed decay has the form:
  %
  %	g2-1 = b + | A * exp ( - K1 * t + K2 / 2 * t*t - K3 / 6 * t*t*t ) |^2
  %
  % NB: here it is |g1|^2 the fitted function
  %
  % The input syntax is simple: just write which order you prefer to use.

  %============================================================================
  % FIT GAMMA vs Q^2
  %============================================================================
  function fit_gamma_q2 ( obj, gammaname, varargin )
  % fit Gamma against Q^2 with a linear (*not* affine) dependence
  % moreover, average A for all angles
  % accepted optional arguments: plot ('yes' or 'no')

   try
    gammaname;
   catch
    error('please select a Gamma you want to fit against q^2 (e.g. Gamma, Gamma1, Gamma2)');
   end

   % find the kind of gamma and label the diffusion constant accordingly
   switch gammaname

    case 'Gamma'
     prop_suff = '';
    case {'Gamma1','Gamma2','Gammae','Gammas'}
     prop_suff = gammaname(6);
    case {'K1','k1'}
     prop_suff = '_K1';
    otherwise
     error('Gamma not recognized!');

   end

   % give some names to the output parameters
   Dname	= ['D',prop_suff];
   dDname	= ['dD',prop_suff];
   Unit_Dname	= ['Unit_D',prop_suff];
   Aname	= ['A',prop_suff];
   A_avname	= ['A',prop_suff,'m'];
   dA_avname	= ['dA',prop_suff,'m'];

   % interpret optional arguments
   options = struct(varargin{:});

   % plot fit results if not forbidden by the user
   try
    options.Plot;
   catch
    options = setfield(options,'Plot','no');
   end

   % initialize some fit options
   fit_function	= '1e6 .* D .* q2';	% unit conversion from ns^-1 to ms^-1
   lowerbond		= 1e-6;
   upperbond		= 1e3;
   startpoint		= 7;

   % fit every concentration
   for i = 1 : length(obj.C)
 
    % this is subtle: I'm using elementwise logical operators and dynamic structure fields
    % In the end we get an array of 1 or 0 indicating us which elements of, say, Gamma1
    % are to be fitted now. I exclude NaN values (bad fits of the correlogramms).
    inde = ~( obj.Data.C - obj.C(i) ) & ~isnan(obj.(gammaname)) & ~isnan(obj.(['d',gammaname]));

    % prepare the vectors for the fit
    q = obj.Data.Q(inde);
    gamma = obj.(gammaname)(inde);
    dgamma = obj.(['d',gammaname])(inde);

    % take also the amplitudes in case of double fit, since they say something
    a = obj.(Aname)(inde);
    da = obj.(['d',Aname])(inde);
 
    % set fit options
    weights =  1 ./ ( dgamma .^ 2 );
    fit_options = fitoptions('Method','NonLinearLeastSquares',	...
       		'Weights', weights, 'MaxFunEvals',100000,	...
           		'Lower',lowerbond,			...
       		'Upper',upperbond,				...
       		'Startpoint',startpoint					);
    fit_type = fittype(fit_function,				...
       		'independent','q2', 				...
       		'coefficients',{'D'},				...
       		'options',fit_options					);

    % perform the numerical
    [ cf gof ] = fit( q' .^2, gamma', fit_type );
     
    % get fit result
    par = coeffvalues (cf); 
    confidence = confint(cf,0.95);
    dpar =  0.5*(confidence(2,:)-confidence(1,:));  

    % output coefficients
    D(i)		= par(1);
    dD(i)		= dpar(1);

    % average amplitudes
    A(i)		= mean(a);
    dA(i)		= sqrt(mean(da .^2));		% assume they have the same mean etc.

   end
 
   % add the D, dD and Unit_D prop to the class
   obj.check_add_prop(Dname,D,dDname,dD,Unit_Dname,'A^2/ns',A_avname,A,dA_avname,dA);

   % print the diffusion constants
   fprintf(['The diffusion constants are: ', Dname, ' ', mat2str(obj.(Dname),2),' ',obj.(Unit_Dname),'\n']);

   % plot fit results if the option is positive
   if strcmp(options.Plot,'yes')
    obj.plot_gamma_q2(gammaname);
   end

  end	% fit_gamma_q2

  %============================================================================
  % CALCULATE HYDRODYNAMIC RATIO
  %============================================================================
  function calculate_H ( self, Dname, sls )
  % this function calculates the H function from SLS and DLS results, using the following formula:
  %
  %	H := D1 * S(q) / D0,	in this version, H:= D1 / ( Kc/R * LIT.M * LIT.D0 )
  %
  % Input parameters:
  % - Dname:	the name of the diffusion constant you want to plot
  % - sls:	this is the SLS class from which one takes the Kc/R values.

   % check that both the diffusion constant exists and sls is a valid SLS class
   try assert( strcmp(class(sls),'SLS') && ismember(Dname,properties(self)) )
   catch error('Your diffusion constant does not exist, or your SLS class is not valid.');
   end
 
   % calculate H ( long vector )
   for i = 1 : length(self.C)
    index_dls	= ( self.Data.C == self.C(i) );				% calculate the dls index
    index_sls	= ( sls.Data.C == self.C(i) );				% calculate the sls index

    D		= mean( self.(Dname)		(index_dls) );		% calculate D TODO: INCORRECT!
    dD		= mean( self.(['d',Dname])	(index_dls) );		% calculate dD TODO: INCORRECT!
    KcR		= mean( sls.Data.KcR		(index_sls) );		% calculate KcR TODO: INCORRECT!
    dKcR	= mean( sls.Data.dKcR		(index_sls) );		% calculate dKCR TODO: INCORRECT!

    H(i)	= D / ( KcR * LIT.BSA.M * LIT.BSA.D0 );			% calculate H
    dH(i)	= H(i) * sqrt( ( dD / D )^2 + ( dKcR / KcR )^2 );	% Gaussian error propagation
   end

   self.check_add_prop(	'H',	H,	...
			'dH',	dH,	...
			'Unit_H','no units');				% add H to the class
   end	% plot_H

  %============================================================================
  % CORRELATE GAMMA
  %============================================================================
  function correlate_gamma ( self, gamma1, gamma2, varargin )
  % this function plots gamma1 vs gamma2, to look for correlations
  % moreover, it calculates the cross-correlation and shows it
  % Accepted optional arguments:
  % - Color:	Base color for nuances, such as Green, Red, etc.
  % - Figure:	Use existing figure or create new one
  % - Type:	type of plot, such as Fill, Points, Lines

   if nargin == 1							% fallback gammas: Gamma1 and Gamma2
    gamma1 = 'Gamma1';	gamma2 = 'Gamma2';
   end

   try		self.(gamma1);	self.(gamma2);				% check for existence of the gammas
   catch	error('Something is wrong with your gammas.');
   end

   options	= struct(varargin{:});					% get optional parameters

   try color = options.Color;	catch	color = 'Green'; end		% color...
   try type = options.Type;	catch	type = 'Fill'; end		% ...type of plot...

   try parameter = options.Parameter; catch parameter = 'C'; end	% ...get the parameter for color nuances...
   nuances = self.get_color(color,length(self.(parameter)));		% ...nuances for plot

   D1_a		= 1e-6 * self.(gamma1) ./ ( self.Data.Q.^2 );		% get the data from the class
   D2_a		= 1e-6 * self.(gamma2) ./ ( self.Data.Q.^2 );
   dD1_a	= 1e-6 * self.(['d',gamma1]) ./ ( self.Data.Q.^2 );
   dD2_a	= 1e-6 * self.(['d',gamma2]) ./ ( self.Data.Q.^2 );

   try fig = options.Figure;
   catch
    fig = self.create_figure;						% create the figure
    xlabel('Faster diffusion [ A^2 / ns ]');				% x label
    ylabel('Slower diffusion [ A^2 / ns ]');				% y label
   
   end
   ax = get(fig,'CurrentAxes');

   for i = 1 : length(self.(parameter))					% group the points for colors, convhull etc.
    D1{i}	= D1_a( self.Data.(parameter)	== self.(parameter)(i) );
    D2{i}	= D2_a(	self.Data.(parameter)	== self.(parameter)(i) );
    dD1{i}	= dD1_a(self.Data.(parameter)	== self.(parameter)(i) );
    dD2{i}	= dD2_a(self.Data.(parameter)	== self.(parameter)(i) );

    switch type								% plot differently according to type
     case 'Points'							% only the fit points
      plot(ax, D1{i}, D2{i},	'*',				...
			'Color',	nuances(i,:),		...
			'MarkerSize',	0.6*self.MarkerSize	);

     case 'Lines'							% only the convex hull
      K	= convhull(D1{i},D2{i});					% calculate the convex hull
      plot(ax, D1{i}(K), D2{i}(K),	'-',			...
			'Color',	nuances(i,:),		...
			'LineWidth',	self.LineWidth		);

     case 'Fill'							% only the area inside the convex hull
      K	= convhull(D1{i},D2{i});					% calculate the convex hull
      p = patch( D1{i}(K), D2{i}(K), nuances(i,:) );			% plot filled area of convex hull	
      set(p,'FaceAlpha',.7,'EdgeAlpha',.5);				% set transparency

     otherwise
      error('Plot type not understood!');
    end

   end

   CC	= ( mean(D1_a .* D2_a) / ( mean(D1_a) * mean(D2_a) )  ) - 1;

   text('Units','normalized','Position',[0.67 0.05],		...
	'FontSize',	self.FontSize,				... 
	'FontName',	self.FontName,				...
	'FontWeight',	self.FontWeight,				...
	'String',['g_1(\',gamma1,',\',gamma2,')',	...
		'=',num2str(CC,3)]);					% show cross-correlation

  end	% correlate gamma

  %============================================================================
  % LOG LIN TOGGLE
  %============================================================================
  function loglin_toggle ( obj, angle, axes)
  % convert a lin_log plot into a log_lin plot, with reasonable limits for axes

   if nargin == 1
    axes = gca;
    angle = obj.Angles(1);
   elseif nargin == 2
    axes = gca;
   end

   % behaviour depends on the actual situation
   if strncmp(get(axes,'xscale'),'log',3)
    set(axes,'yscale','log','xscale','lin');
    % xmax depends strongly on the angle
    xmax	= 0.11 + 0.6 * ( obj.Angles(length(obj.Angles)) - angle ) .^2 / ( obj.Angles(length(obj.Angles)) - obj.Angles(1) ) .^2;
    axis(axes, [ 0 xmax 1e-3 1 ] );
   else
    set(axes,'yscale','lin','xscale','log');
    axis(axes,'auto' );
   end

  end	% end of function

 end	% end of public methods
 
 %============================================================================
 % PRIVATE METHODS
 %============================================================================
 methods ( Access=private )

 end	% end of private methods
end	% end of class
