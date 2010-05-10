%============================================================================
%
% DLS_gen class for grouping various information from the single DLS classes
% and present them coherently
%
% Input is an array of DLS classes
%
%============================================================================

classdef LIT
 properties ( Access = public, Constant )					% with the 'Constant' attribute, they can be called from outside without an instance

 % conversion factors
 Na		=	6.02214179e23;						% Avogadro number
 u_in_Kg	=	1.660538782e-27;					% unit of atomic mass
 kb		=	1.3806504e-23;						% Boltzmann constant in J/K
 T0		=	-273.15;						% absolute zero in °C

 % BSA properties
 BSA = struct(		'M',			66411,			...
			'Unit_M',		'Da',			...
			'dndc_H2O',		[ 183 181 ],		...
			'dndc_100mM_NaCl',	180,			...
			'Unit_dndc',		'mul/g',		...
			'dndc_lambda',		'632 nm'			);

 % buffers properties
 CH3COOH_M_g_mol	=	60.05;
 NaCH3COO_M_g_mol	=	82.03;



 %end of properties
 end

 methods

  % constructor
  function lit = LIT ()

   fprintf('\nWelcome! The syntax is the usual: <LIT>.<property>.\n\n');
 
  % end of constructor
  end

 % end of methods
 end

% end of LIT class
end
