function plot ( varargin )

 a	= varargin;				% options for the fit

 if ishandle(a{1}) & ~ismember(class(a{1}),{	'SLS.Sample',		'Sample',	...
						'SLS.Experiment',	'Experiment'	})
  ax = a{1};	a(1) = [];
 else
  ax = axes;
 end
 hold(ax,'on');

 obj	=	a{1};

 for i = 1 : length(obj)			% obj could be a class array

  point	=	obj(i);
 
  x	=	[ point.(a{2}) ];
  y	=	[ point.(a{3}) ];
 
  [ x i ]	= sort(x);
  y		= y(i);
 
  try	s = a{4};
  catch	s = 'o';
  end
 
  options = a(5:length(a));

  plot(ax,	x,	y,	s,	...
 		'LineWidth',	4,	...
 		'MarkerSize',	10,	...
 		options{:}		);
 end

 legend('show');

end
