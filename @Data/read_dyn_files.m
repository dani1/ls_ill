function [ Angles Tau G dG ] =  read_dyn_files( path, runstart, runs )
% read text files to get DLS data

 % read every file
 for i = 1 : runs

  % define the file names
  dyn_files{i}	= [path,'00',num2str(runstart+i-1,'%2.2u'),'.ASC'];

  angle = [];
  tau	= [];
  g	= [];
  dg	= [];
 
  fid = fopen ( dyn_files{i} );		% open the dynamic file
 
  % read the angle
  str = fgetl(fid);
  while ~strcmp(str,'"Correlation"')
   if strfind(str,'Angle')
    [ tmp tmp tmp angle ] = strread(str, '%s %s %s %f'); 
   end
   str = fgetl(fid);
  end
 
  % read tau and g
  while ~strcmp(str,'')
   str=fgetl(fid);
   [ t_t g_t tmp tmp tmp ] = strread(str, '%f %f %f %f %f');
   tau		= [ tau	; t_t ];
   g		= [ g	; g_t ];
  end
 
  % read stddev
  while ~strcmp(str,'"StandardDeviation"') str=fgetl(fid); end
  while ~feof(fid)
   str=fgetl(fid);
   [ tmp dg_t ]	= strread(str, '%f %f');
   dg		= [ dg	; dg_t ];
  end
 
  fclose(fid);				% close the dynamic file

  % store the info
  Angles(i)	= angle;
  Tau{i}	= tau;
  G{i}		= g;
  dG{i}		= dg;

 end

end	% read_dyn_files
