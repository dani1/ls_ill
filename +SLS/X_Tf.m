function [ X_T f ] = X_Tf( T, C, KcR )
% the osmotic isothermal compressibility is related to the structure factor at
% vanishing Q by the following formula:
%
%	X :=	1/c * ( dc / dP )_T	= 1 / [ Na * kT * c * M * S(c, q -> 0)	]
%					= 1 / [ Na * kT * c * 	(Kc/R) 		]
%					= 1 / [ RT * c * 	(Kc/R) 		]
%
% where X is the compressibility, c is the concentration in mass/volume, N_A is
% Avogadro's number and S(c,q) is the static structure factor.
%
% Please note that X(c) is independent on Q (because it is defined for Q-->0) and on M.

 Na	= Constants.Na;
 kb	= Constants.kb;

 f	= inline(	'1 ./ ( Na .* kb .* T .* C .* KcR )', 'Na','kb','T','C','KcR');

 X_T	= f(Na,kb,T,C,KcR);

end
