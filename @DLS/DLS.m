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
   dls.Data.Tau		= horzcat(dataclasses.Tau);
   dls.Data.G		= horzcat(dataclasses.G);
   dls.Data.dG		= horzcat(dataclasses.dG);

   % create unique properties with the update_uniques function in the Utils class
   dls.update_uniques('C','Q','Angles');

   % monitor properties with listeners in the Utils class
   dls.monitorprop('C','Q','Angles');

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

   try
    method;
   catch
    method = 'Double';
   end
 
   % interpret the optional arguments
   options = struct(varargin{:});

   % perform the gamma vs q^2 fit if not prohibited by the user
   try
    options.fit_gamma_q2;
   catch
    options = setfield(options,'Fit_gamma_q2','yes');
   end

   % set the fit parameters, according to the fit method chosen
   % general parameters (background)
   fit_common	= 'BKG';
   coeff_common	= {'BKG'};
   lower_common	= 0;
   upper_common	= 1e-1;
   start_common	= 1e-3;
 
   switch method
 
    case 'Single'
     fit_function	= [ fit_common,'+ A * exp( - 2 * Gamma * t )'	];
     coefficients	= [ coeff_common	{'A'}	{'Gamma'}	];
     lowerbond		= [ lower_common	0.8	1e-2		];
     upperbond		= [ upper_common	1.2	1e5		];
     startpoint		= [ start_common	1.1	1e1		];
 
    case 'Double'
     fit_function	= [ fit_common,'+ ( A1 * exp( - Gamma1 * t ) + A2 * exp( - Gamma2 * t ) ).^2'	];
     coefficients	= [ coeff_common	{'A1'}	{'Gamma1'}	{'A2'}	{'Gamma2'}		];
     lowerbond		= [ lower_common	0	1e-1		0	1e-3			];
     upperbond		= [ upper_common	1	1e9		1	1e2			];
     startpoint		= [ start_common	0.3	1e2		0.3	1e-1			];
 
    case 'Single+Streched'
     fit_function	= [ fit_common,'+ ( Ae * exp( - Gammae * t ) + As * exp(-( Gammas *t).^b )) .^2'];
     coefficients	= [ coeff_common	{'Ae'}	{'Gammae'}	{'As'}	{'Gammas'}	{'b'}	];
     lowerbond		= [ lower_common	0	1e-1		0	1e-3		0	];
     upperbond		= [ upper_common	1	1e9		1	1e2		1	];
     startpoint		= [ start_common	0.3	1e2		0.3	1e-1		0.8	];
 
    otherwise
     error('Method not recognized!');
   end
 
   % fit for all data
   for l = 1 : length(obj.Data.Tau)
 
     % select which range of data to use for the fit. The problem is mainly
     % that for long lag times I'm fitting background, whereas for very short
     % lag times it is afterpulse noise of the APDs
% TO DO: control the behaviour of the correlations at 100 times higher than the
%        first decay time. If there is still something important, do not fit
%        the data and impose NaN for both Gamma and dGamma. One can not live with
%        scorpions is his bed!
     inde	= obj.Data.Tau{l} > 50e-6 & 1;
     tau	= obj.Data.Tau{l}(inde);
     g		= obj.Data.G{l}(inde);
     dg		= obj.Data.dG{l}(inde);


     % set the weights for the fit
     weights =  1 ./ dg .^ 2;

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

   % add the parameters to the class
   for k = 1 : length(coefficients)
    obj.check_add_prop(coefficients{k},coout(k,:));
    obj.check_add_prop(['d',coefficients{k}],dcoout(k,:));
   end
   obj.check_add_prop('Unit_Gamma',[ obj.Unit_Tau,'^{-1}']);

   % if the option fit_gamma_q2 is 'yes' (default), fit linearly gamma to Q^2
   if strcmp(options.Fit_gamma_q2,'yes');

     % fit all gammas
     for k = 1 : length(coefficients)
      if regexp(coefficients{k},'Gamma')
       obj.fit_gamma_q2(coefficients{k});
      end
     end

   end

  end	% fit_correlations
  
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
   obj.check_add_prop(Dname,D,dDname,dD,'Unit_D','A^2/ns',A_avname,A,dA_avname,dA);

   % print the diffusion constants
   fprintf(['The diffusion constants are: ', Dname, ' ', mat2str(obj.(Dname),2),' ',obj.Unit_D,'\n']);

   % plot fit results if the option is positive
   if strcmp(options.Plot,'yes')
    obj.plot_gamma_q2(gammaname);
   end

  end	% fit_gamma_q2

  %============================================================================
  % PLOT GAMMA vs Q^2
  %============================================================================
  function plot_gamma_q2 ( obj, gammaname, varargin )
  % plot the Gammas from the fit against Q^2 and, if found, the corresponding diffusion function
  % accepted optional arguments: C, 

   try
    gammaname;
   catch
    error('please select a Gamma you want to fit against q^2 (e.g. Gamma, Gamma1, Gamma2)');
   end

   options	= struct(varargin{:});						% eat the optional arguments as a struct

   try										% concentrations
    options.C;
   catch
    options	= setfield(options,'C',obj.C);
   end

   for i = 1 : length(options.C)

    % this is subtle: I'm using elementwise logical operators and dynamic structure fields
    % In the end we get an array of 1 or 0 indicating us which elements of, say, Gamma1
    % are to be fitted now. I exclude NaN values (bad fits of the correlogramms).
    inde = obj.Data.C == options.C(i) & ~isnan(obj.(gammaname)) & ~isnan(obj.(['d',gammaname]));

    % prepare the vectors for the fit
    q = obj.Data.Q(inde);
    gamma = obj.(gammaname)(inde);
    dgamma = obj.(['d',gammaname])(inde);

    obj.figure_gamma_q2;							% create the figure
 
    % title
    titlename = ['Plot of ', char(gammaname),' from the fit at cp = ', num2str(options.C(i)),' againt q^2'];
    title(titlename,'FontSize',14);

    % plot with error bars
    errorbar (q .^2 , gamma, dgamma, 'o','MarkerSize', 10 );

    % check if the right diffusion constant has been fitted, and if yes plot it
    try
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
     D		= obj.(Dname)(obj.C == options.C(i));
     % this is the actual try
     Dplo = plot(q .^2, 1e6 * q .^2 * D, '-', 'LineWidth', 2 );
     
     legend(Dplo,[Dname,' = ',num2str(D,2),' ',obj.Unit_D],'Location','NorthWest');
    end

   end

  end	% plot_gamma_q2

  %============================================================================
  % PLOT D vs C
  %============================================================================
  function plot_D ( obj, Dname, varargin )
  % accepted optional arguments: Figure

   try
    Dname;
   catch
    error('Please select which diffusion constant to plot.');
   end

   % interpret optional arguments
   options = struct(varargin{:});

   % if a figure is specified, use it instead of a new one
   try
    fig = options.Figure;
   catch
    fig = obj.figure_D_C;
%    title([Dname,' versus protein concentration'],'FontSize',14);
   end

   % plot with errorbars
   plot_options = { '-o','LineWidth',4,'MarkerSize',10 };

   % add color to plot_options if chosen by the user
   try
    plot_options = [ plot_options, 'Color', options.Color ];
   end

   errorbar(get(fig,'CurrentAxes'),obj.C,obj.(Dname),obj.(['d',Dname]), plot_options{:});

  end	% plot_D

  %============================================================================
  % PLOT H
  %============================================================================
  function plot_H ( obj, Dname, sls )
  % this function plots the H function (in this version, only something proportional to it)
  % from SLS and DLS results, using the following formula:
  %
  %	H := D1 * S(q) / D0	( in this version, H:= D1 / Kc/R * M_lit / D1(1) )
  %
  % Input parameters:
  % - Dname:	the name of the diffusion constant you want to plot
  % - sls:	this is the SLS class from which one takes the Kc/R values.

  % check that both the diffusion constant exists and sls is a valid SLS class
  if ( strcmp(class(sls),'SLS') && ismember(Dname,properties(obj)) )

   % mean over the Q values for the SLS results
   for i = 1 : length(obj.C)
    inde	= sls.Data.C==obj.C(i);
    H(i)	= obj.(Dname)(i) / ( sum( sls.Data.KcR(inde) ./ ( sls.Data.dKcR(inde) .^2 ) )	/ sum ( 1 ./ sls.Data.dKcR(inde) .^2 ) ) / 67000 / obj.(Dname)(1);
    dH(i)	= obj.(Dname)(i) / ( sum( sls.Data.dKcR(inde) ./ ( sls.Data.dKcR(inde) .^2 ) )	/ sum ( 1 ./ sls.Data.dKcR(inde) .^2 ) ) / 67000 / obj.(Dname)(1);
   end
 
   fig = figure;
   hold all;
   set( gcf, 'color', 'white' );
   set( gca,	'box', 'on',						...
		'FontSize',36,'FontName','Courier','FontWeight','bold',	...
		'XMinorTick','on', 'YMinorTick','on');
   xlabel(['Protein concentration [ ',obj.Unit_C,' ]'],'FontSize',36);
   ylabel(['H [ a.u. ]'],'FontSize',36);

   plot(obj.C,H,'o-','LineWidth',4);

  end

  end	% plot_H

  %============================================================================
  % CORRELATE GAMMA
  %============================================================================
  function correlate_gamma ( self, gamma1, gamma2, varargin )

   options	= struct(varargin{:});					% get optional parameters

   try color = options.Color;	catch	color = 'Green'; end		% color
   nuance = self.get_color(color,1);					% nuance for plot

   try		self.(gamma1);	self.(gamma2);				% check for existence of the gammas
   catch	error('Something is wrong with your gammas.');
   end

   D1_a		= 1e-6 * self.(gamma1) ./ ( self.Data.Q.^2 );			% get the data from the class
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

   plot(ax, D1_a, D2_a,	'*',					...
			'Color',	nuance,			...
			'MarkerSize',	0.6*self.MarkerSize	);	% plot

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

  %============================================================================
  % create a figure for gamma vs q2
  %============================================================================
  function fig = figure_gamma_q2 ( obj )
   
   fig = figure;
   hold all;
   set( gcf, 'color', 'white' );
   set( gca,	'box', 'on',						...
		'FontName','Courier','FontSize',36,'FontWeight','bold',	...
		'XMinorTick','on', 'YMinorTick','on');
   xlabel(['Q^2 [ ',obj.Unit_Q,'^2 ]'],'FontSize',36);
   ylabel(['\Gamma [ ',obj.Unit_Tau,'^{-1} ]'],'FontSize',36);
 
  end	% end of figure function

  %============================================================================
  % create a figure for D vs C
  %============================================================================
  function fig = figure_D_C ( obj )
   
   fig = figure;
   hold all;
   set( gcf, 'color', 'white' );
   set( gca,	'box', 'on',						...
		'FontSize',36,'FontName','Courier','FontWeight','bold',	...
		'XMinorTick','on', 'YMinorTick','on');
   xlabel(['Protein concentration [ ',obj.Unit_C,' ]'],'FontSize',36);
   ylabel(['Diffusion constant [ ',obj.Unit_D,' ]'],'FontSize',36);
 
  end	% end of figure function

 end	% end of private methods
end	% end of class
