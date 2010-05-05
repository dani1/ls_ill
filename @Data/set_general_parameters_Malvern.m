%=========================================================================================
% low-level routine for giving general data about the Malvern instrument
% Malvern is so user-friendly that the instrument does not say anything
% Thus, we just generate the information according to our knowledge
%=========================================================================================
function [ lambda unit_lambda n_set ] = set_general_parameters_Malvern

 lambda		= 632;				% is it a HeNe laser?
 unit_lambda	= 'nm';

 n_set		= 1.332;			% assuming a water standard?

end	% read_general_file_ALV
