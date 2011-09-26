%============================================================================
%
% BSA class
% 
% This class contains information about the BSA molecule. I can present
% results coherently using these values.
%
%============================================================================
classdef BSA
 properties ( Access = public, Constant )					% with the 'Constant' attribute, they can be called from outside without an instance

 M			=	66411
 Unit_M			=	'Da'
 dndc_H2O		=	[ 183 181 ]
 dndc_100mM_NaCl	=	180
 Unit_dndc		=	'mul/g'
 dndc_lambda		=	'632 nm'
 D0			=	6.0
 Unit_D0		=	'A^2/ns'
 v0			=	1.140e-3
 Unit_v0		=	'l/g'
 vmol			=	155
 Unit_vmol		=	'nm^3'
 R			=	3.60
 Unit_R			=	'nm'

 % buffers properties
 CH3COOH_M_g_mol	=	60.05;
 NaCH3COO_M_g_mol	=	82.03;

 %end of properties
 end

% end of LIT class
end
