%=========================================================================================
% low-level routine for reading general data from autocorrelation files
% the first file is opened and info about instrument and lambda are retrieved
%=========================================================================================
function [ lambda unit_lambda n_set ] = read_general_file_ALV ( genpath, runstart )

 if nargin == 1
  runstart	= 0;
 end

 genpath	= [ genpath num2str(runstart,'%4.4d') '.ASC' ];
 
 fid = fopen( genpath );								% open the file
 
 while ~feof(fid)									% till the end of file, do the following

  str = fgetl(fid);									% read every line
  if strfind(str,'Refractive') [ tmp tmp n_set ] = strread(str, '%s %s %f');	end	% n_set
  if strfind(str,'Wavelength')								% wavelength
   [ tmp ul tmp lambda ] = strread(str, '%s %s %s %f');
   % try to interpret the unit of lambda
   ul=char(ul);
   if strfind(ul,'nm')		unit_lambda = 'nm';
   elseif strfind(ul,'A')	unit_lambda = 'A';				end
  end

 end
 
 fclose(fid);										% close the file

end	% read_general_file_ALV
