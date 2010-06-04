%============================================================================
%
% Graphics class for generating and controlling graphic objects
% 
% The idea is to collect some functions useful across all LS classes,
% which become subclasses of this one and inherit the methods
%
%============================================================================

classdef Graphics < handle
 %===========================================================================
 % PUBLIC PROPERTIES
 %===========================================================================
 properties ( Access = public, Hidden )						% subclasses inherit these

  FontSize	= 36
  FontName	= 'Courier'
  FontWeight	= 'bold'
  XMinorTick	= 'on'
  YMinorTick	= 'on'

  MarkerSize	= 10
  LineWidth	= 4

 end	% properties

 %===========================================================================
 % METHODS
 %===========================================================================
 methods ( Access = public )							% subclasses inherit the methods

  %==========================================================================
  % CREATE FIGURE
  %==========================================================================
  function [ fig ax ] = create_figure ( self )
  % This function creates a figure without axes
   
   fig	= figure;								% create the figure
   ax	= gca;									% get the axes
   hold all;									% hold the axes for plots
   set( gcf, 'color', 'white' );						% background
   set( gca,	'box', 'on', 				...
		'FontSize',	self.FontSize,		...
		'FontName',	self.FontName,		...
		'FontWeight',	self.FontWeight,	...
		'XMinorTick',	self.XMinorTick,	...
		'YMinorTick',	self.YMinorTick		);

  end	% create_figure

  %==========================================================================
  % GET COLOR FOR PLOTS
  %==========================================================================
  function nuances = get_color ( self, color, len )
  % function to give the nuances of colors for the plots as a function of the
  % color itself and the number of nuances

   R = zeros(1,len);
   G = zeros(1,len);
   B = zeros(1,len);
 
   switch color
    case 'Green'			% Green
     R(:) = 0.1;
     G = [1:len]/(len+1);
     B(:) = 0.2;
 
    case 'Red'				% Red/Yellow
     R(:) = 0.8;
     G = [1:len]/(len+1);
     B(:) = 0.2;
 
    case 'Blue'				% Blue
     R(:) = 0.1;
     G = [1:len]/(len+1);
     B(:) = 0.8;
 
    case 'Pink'				% Pink/Lime
     R(:) = 0.7;
     G = [1:len]/(len+1);
     B(:) = 0.6;
 
    case 'Gray'				% Grayscale
     R = [1:len]/(len+1);
     G = [1:len]/(len+1);
     B = [1:len]/(len+1);

    otherwise
     error('Color not recognized!');
   end
 
   nuances	= [ R' G' B' ];

  end	% get_color

  %============================================================================
  % SET XLABEL ACCORDING TO INDEPENDENT
  %============================================================================
  function xl = set_xlabel ( self, independent )

   switch independent
    case 'C'							% concentration
     xl = ['Protein concentration [ ',self.Unit_C,' ]'];

    case 'Phi'							% volume fraction
     xl	= ['\Phi'];

    case 'I'							% ionic strength
     xl = ['ionic strength [ ',self.Unit_I,' ]'];

    case 'Q2'							% Q squared
     xl = ['Q^2 [ ',self.Unit_Q2,' ]'];

    otherwise							% not recognized?!
     error('Independent not recognized!');
   end

  end	% set_xlabel

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
   [ax, x, y, dy, color] = self.prepare_errorbar( name, optargs{:} );			% prepare the plot
   self.make_errorbar ( ax, x, y, dy, color );						% make the plot

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
   [ ax, x, y, color ] = self.prepare_errorbar ( name, optargs{:} );			% prepare the plot
   self.make_plot ( ax, x, y, color );							% make the plot

  end % plot

  %==========================================================================
  % PREPARE PLOT
  %==========================================================================
  function [ ax, x, y, color ] = prepare_plot ( self, name, varargin )
  % this function prepares the plot for the overloaded plot method
  % note that 'plot' could be in turn overloaded by subclasses, to handle
  % class-specific properties such as compressibility, etc.

   options	= struct(varargin{:});							% get optional arguments

   try		y	= self.Data.(name);						% try looking in the Data struct
   catch try	y	= self.(name);							% otherwise, look directly
   catch	error(['Property ',name,' not found!']); end				% print an error if the prop is not found
   end

   try	color = options.Color;			catch color = 'Green';		end	% default color: Green

   try	independent = options.Independent;	catch independent = 'C';	end	% ...default independent: C

   if length(y) == length(self.Data.C)	x	= self.Data.(independent);		% if the property is long, use the Data independent...
   else					x	= self.(independent);		end	% ...else, use the unique independent

   try											% try using an existing figure...
    fig		= options.Figure;
    ax		= get(fig,'CurrentAxes');
   catch										% ...or create a new one
    [fig ax]	= self.create_figure;

    try [ xunit yunit ] = self.get_units( independent, name );			end	% get the units for the axes labels

    try											% try to set an x label
     if isempty(xunit)	xlabel( ax,	independent);
     else		xlabel( ax,	[independent,' [ ',xunit,' ]']);	end
    end

    try											% try to set an x label
     if isempty(yunit)	ylabel( ax,	name);
     else		ylabel( ax,	[name,' [ ',yunit,' ]']);		end
    end

   end

  end	% prepare_plot


  %==========================================================================
  % PREPARE ERRORBAR
  %==========================================================================
  function [ ax, x, y, dy, color ] = prepare_errorbar ( self, name, varargin )

   optargs	= varargin;

   [ ax, x, y, color ] = self.prepare_plot( name, optargs{:} );

   try		dy	= self.Data.(['d',name]);					% try looking in the Data struct
   catch try	dy	= self.(['d',name]);						% otherwise, look directly
   catch	error(['Property ',['d',name],' not found!']); end			% print an error if the prop is not found
   end

  end	% prepare_errorbar

  %==========================================================================
  % GET AXES UNITS
  %==========================================================================
  function [ xunit yunit ] = get_units( self, independent, dependent )
  % this function eats the axes names and outputs the units for labels
  % NB: this method can be overloaded by subclasses

   xunit	= self.(['Unit_',independent]);					% try to get the units of x axis...
   yunit	= self.(['Unit_',dependent]);					% try to set the units for y axis...

  end	% get_units

  %==========================================================================
  % MAKE PLOT
  %==========================================================================
  function make_plot ( self, ax, x, y, color )
  % this function makes the plot after that it has been prepared

   if length(y) < length(self.Data.C)							% if the property is short...

    nuance	= self.get_color(color,1);						% ...get the nuance for the plot...
    plot( ax, x, y, 'o-',			...
		'Color',	nuance,			...
		'MarkerSize',	self.MarkerSize,	...
		'LineWidth',	self.LineWidth		);				% ...and plot!

   else											% ...otherwise, use the unique independent

    try	average = options.Average; 		catch average = 'yes'; end		% default: average
    switch	average									% act differently if one averages or not

     case 'yes'										% in this case, do only one plot

      [ xm ym ] = self.average( x, y );							% get the average values

      nuance	= self.get_color(color,1);						% get the nuance for the plot
      plot( ax, xm, ym, 'o-',			...
		'Color',	nuance,			...
		'MarkerSize',	self.MarkerSize,	...
		'LineWidth',	self.LineWidth		);				% plot

     case 'no'										% in this case, color the lines differently parametrically

      try parameter = options.Parameter;						% try to get the parameter...
      catch parameter = self.set_parameter(independent);				% ...or fall back to default
      end

      p	= self.(parameter);								% get the parametric vector

      nuances	= self.get_color(color,length(p));					% get the nuance for the plot
      for i = 1 : length(p)
       index	= ( self.Data.(parameter) == p(i) );					% get the index for this plot
       xp	= x(index);
       yp	= y(index);

       plot( ax, xp, yp, 'o-',			...
		'Color',	nuances(i,:),		...
		'MarkerSize',	self.MarkerSize,	...
		'LineWidth',	self.LineWidth		);				% plot
      end

    end
   end

  end	% make_plot

  %==========================================================================
  % MAKE ERRORBAR
  %==========================================================================
  function make_errorbar ( self, ax, x, y, dy, color )
  % this function makes the errorbar plot after that it has been prepared

   if length(y) < length(self.Data.C)							% if the property is short...

    nuance	= self.get_color(color,1);						% get the nuance for the plot
    errorbar( ax, x, y, dy, 'o-',			...
		'Color',	nuance,			...
		'MarkerSize',	self.MarkerSize,	...
		'LineWidth',	self.LineWidth		);				% plot

   else											% ...else, use the unique independent

    try	average = options.Average; 		catch average = 'yes'; end		% default: average
    switch	average									% act differently if one averages or not

     case 'yes'										% in this case, do only one plot

      [ xm ym dym ] = self.average( x, y, dy );						% get the average values

      nuance	= self.get_color(color,1);						% get the nuance for the plot
      errorbar( ax, xm, ym, dym, 'o-',			...
		'Color',	nuance,			...
		'MarkerSize',	self.MarkerSize,	...
		'LineWidth',	self.LineWidth		);				% plot

     case 'no'										% in this case, color the lines differently parametrically

      try parameter = options.Parameter;						% try to get the parameter...
      catch parameter = self.set_parameter(independent);				% ...or fall back to default
      end

      p	= self.(parameter);								% get the parametric vector

      nuances	= self.get_color(color,length(p));					% get the nuance for the plot
      for i = 1 : length(p)
       index	= ( self.Data.(parameter) == p(i) );					% get the index for this plot
       xp	= x(index);
       yp	= y(index);
       dyp	= dy(index);

       errorbar( ax, xp, yp, dyp, 'o-',			...
		'Color',	nuances(i,:),		...
		'MarkerSize',	self.MarkerSize,	...
		'LineWidth',	self.LineWidth		);				% plot
      end

    end
   end

  end	% make_errorbar

 end	% methods

end	% Graphics class
