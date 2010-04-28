%=========================================================================================
% low-level routine for reading general data from autocorrelation files
% the first file is opened and info about instrument and lambda are retrieved
%=========================================================================================
function [ instrument lambda unit_lambda n_set ] = load_general_data( path, runstart )
 
 first_file = [path,'00',num2str(runstart,'%2.2u'),'.ASC'];

 if exist(first_file, 'file') ~= 2
  error('I can not find the correlation files! Please check your path.');
 else
  fid = fopen( first_file );	% open the file
 
  while ~feof(fid)

   str = fgetl(fid);
   if strfind(str,'ALV')	instrument=str;					end	% instrument
   if strfind(str,'Refractive') [ tmp tmp n_set ] = strread(str, '%s %s %f');	end	% n_set
   if strfind(str,'Wavelength')								% wavelength
    [ tmp ul tmp lambda ] = strread(str, '%s %s %s %f');
    % try to interpret the unit of lambda
    ul=char(ul);
    if strfind(ul,'nm')		unit_lambda = 'nm';
    elseif strfind(ul,'A')	unit_lambda = 'A';				end
   end

  end
  
  fclose(fid);			% close the file
 end

end	% load_general_data
