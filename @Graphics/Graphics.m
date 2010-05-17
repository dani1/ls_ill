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

 end	% methods

end	% Graphics class
