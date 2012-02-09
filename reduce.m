function [ xr yr dyr ] = reduce ( x, y, dy )
% average the values of y according to x

 if nargin < 3
  dy = ones(1,length(y));
 end 

 x	= x(:);
 y	= y(:);
 dy	= dy(:);

 w	= 1 ./ dy.^2;
 xr	= unique(x);
 lr	= length(xr);

 f	= @(x,w)(sum(x.*w) / sum(w)); 

 for i = 1 : lr
  ind	= ( x == xr(i) );

  yr(i)		= f( y(ind),	w(ind) );
  dyr(i)	= f( dy(ind),	w(ind) );

 end
 
end
