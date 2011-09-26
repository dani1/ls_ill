% MULTI-LORENTZ EXAMPLE
%x   = -5 : 0.01 : 5; 
%y   = 3 * 1/pi*0.4./(x.^2 + 0.4^2) + 5*1/pi*2./(x.^2 + 2^2); 
%dy  = 0.25 * randn(1, length(y)); 
%var = 0.25^2*ones(1, length(y)); 

% MULTI_EXPONENTIAL EXAMPLE
x	= 1 : 0.5 : 20;
y	= 0.5 * exp( - x ./ 2 ) + 0.5 * exp( -x ./ 4 );
%y	= exp( - x ./ 5 );
dy 	= 0.015 * randn(1, length(y)); 
var	= 0.50^2*ones(1, length(y)); 

[s1, g1, b1] = contin(x, y+dy, var, min(x), max(x), 10*length(x), 0.1, 0);

legend('off');
cla;

ax = gca;
hold all;

norm	= max(g1);
errorbar(ax,x,norm*y,norm*dy,	'LineWidth',	3);
plot(ax,s1,g1,			'LineWidth',	3);

xlim([0 100]);
set(gca,'XScale','log');
set(gca,'XScale','log');

legend('Original function',	'Inverse Laplace Transform: SPG'	);
