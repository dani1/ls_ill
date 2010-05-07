function [ ang kcr dkcr ] = read_slspath_Malvern ( slspath )
% read SLS data from files for the Malvern Zetasizer Nano in Tuebingen

  % initialize to void, if the function outputs void then skip SLS
  ang	= [];
  kcr	= [];
  dkcr	= [];

  % initialize to void some other vectors
  stand	= [];
  solv	= [];
  samp	= [];
 
 if exist(slspath,'file') ~= 2					% check for existance
  slspath = input([slspath, ...
	' not found. Please enter the full path of the static file or let empty for skipping: '],...
			's');
 end

 % check whether the new file exists
 if exist(slspath,'file') == 2

  fid	 = fopen(slspath);						% open the file

  while ~feof(fid)							% go down the whole file and do the following

   str	= fgetl(fid);							% get the line

   if strfind(str,'Scattering standard count')				% if this is a measurement of standard...

    [ tmp tmp tmp stand_tmp ]	= regexp( str, '[\d,]+' );		% this sequence is ridiculous, but we need that number
    stand_tmp			= regexprep( stand_tmp, ',','.');
    stand_tmp			= str2num( stand_tmp{:} );

    stand = [ stand stand_tmp ];					% grow the stand string;

   elseif strfind(str,'Solvent only scattering count')			% if this is a measurement of pure solvent...

    [ tmp tmp tmp solv_tmp ]	= regexp( str, '[\d,]+' );		% this sequence is ridiculous, but we need that number
    solv_tmp			= regexprep( solv_tmp, ',','.');
    solv_tmp			= str2num( solv_tmp{:} );

    solv = [ solv solv_tmp ];						% grow the solv string;

   elseif strfind(str,'Sample scattering count')			% if this is a measurement of sample...

    [ tmp tmp tmp samp_tmp ]	= regexp( str, '[\d,]+' );		% this sequence is ridiculous, but we need that number
    samp_tmp			= regexprep( samp_tmp, ',','.');
    samp_tmp			= str2num( samp_tmp{:} );

    samp = [ samp samp_tmp ];						% grow the sample string;

   end
  end

  fclose(fid);			% close the file

 end

 st	= mean(stand);		dst	= sqrt(sum((stand-st).^2)/length(stand));
 so	= mean(solv);		dso	= sqrt(sum((solv-so).^2)/length(solv));
 sa	= mean(samp);		dsa	= sqrt(sum((samp-sa).^2)/length(samp));

 ang	= 173;								% this is fixed: NISB
 I	= ( sa -so ) / st;						% scattering intensity (normalized)
 dI	= I * sqrt(((dsa^2+dso^2)/(sa -so)^2 )+(dst/st)^2);		% Gauss error propagation

 kcr	= 1 / I;							% INCOMPLETE!
 dkcr	= dI / I * kcr;							% INCOMPLETE!

end	% read_slspath_Malvern
