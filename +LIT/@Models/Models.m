%============================================================================
%
% Models class
% 
% This class contains the various theoretical models for fits, etc.
%
%============================================================================
classdef Models
 properties ( Access = public, Constant )					% with the 'Constant' attribute, they can be called from outside without an instance

  Hydro	= struct(	'Disordered_HS',	'(1-Phi).^6.55',					...
			'Ordered_HS',		['(1 - 1.5.*Phi.^(1/3) + 1.5.*Phi.^(5/3) - Phi.^2)',	...
						'./ (1 + 2/3*Phi.^(5/3))']				);		% models for the hydrodynamic function		

 end	% constant properties

end	% models class
