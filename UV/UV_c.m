%==========================================================================
% UV_c
%
% UV/VIS spectra measure the absorbance (Beer-Lambert-law):
% 			A	=	l * eps * c
% where:
% - eps is the molar extinction coefficient (wavelength dependent);
% - l is path length;
% - c the molar concentration.
%
% Proteins strongly absorb at 279nm.
%
% ---	SYNTAX	---
%
% This function takes as input:
% - the path of the absorption file to use;
% - the pipette colors, quantities and the number of protein pipettes,
%   for every dilution step, e.g.
%
% Syntax example:
% 
% [ c dc ] =  UV_c('../meshdata/091127-28-Marcus/BSA_1mgml_200nm-1.jws.txt', {'blue' 'yellow' 'yellow'}, [ 500 100 100 ], 1 )
%
% The output is the mass concentration as grams of solute over liters of solution
%==========================================================================
function [ c dc ] = UV_c ( pathUV, varargin )

 [ gamma dgamma ] = dilution(varargin{:});				% get the dilution factor

 % cuvette specific parameters
 length		= 1;						% [cm]

 % BSA specific parameters
 eps_BSA	= 4.42881 * 1E4;				% [M^-1 * cm^-1]
 mw_BSA		= 66431.5;					% molecular weight of bsa [g/mol]
 
 absorbance	= read_UV_file ( pathUV );			% read text file
     
 % Marcus checked the hardware specifications of the V630 spectrometer 
 % where they stated to have an error of 0.3 for the transmittance   T [%]
 % using the formula T = 10^-A, I derived the error for A
 A	= absorbance( wavelength == 279 );			% get abs at the right lambda
 dT 	= 0.3;
 dA 	= dT / (100 * log(10) * 10^(-A));
     
 c	= A / (eps_BSA * length * gamma ) * mw_BSA;		% concentration [mg/ml]
 dc	= c * sqrt( ( dA / A ) ^2 + ( dgamma / gamma ) ^2 );	% conc error	[mg/ml]
     
 %==========================================================================
 % NESTED FUNCTION: read the data from files
 %==========================================================================
 function absorbance = read_UV_file ( pathUV )
  str = '';  wavelength = []; absorbance = [];			% initialize
  fid = fopen( pathUV );					% open file

  while ~strcmpi(str,'XYDATA')
   str = fgetl(fid);
  end
  while ~strcmp(str,'') && ~feof(fid)
   str	= fgetl(fid);
   [ wl ab ]	= strread(str,'%f %f');				% get lambda and abs
   wavelength	= [ wavelength; wl ];
   absorbance	= [ absorbance; ab ];
  end

  fclose(fid);							% close file
 end	% read_UV_file

end	% UV_c
