function point	= read_dynamic_file ( path )

 fid	= fopen( path );
% disp(path);
 tau	= [];
 g	= [];
 dg	= [];
 
 % read the angle
 str = fgetl(fid);
 while ~strcmp(str,'"Correlation"')
  if strfind(str,'Temperature')
   [ tmp tmp tmp T ] = strread(str, '%s %s %s %f'); 
  elseif strfind(str,'Angle')
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
 while (~strcmp(str,'"StandardDeviation"') && ~feof(fid)) str=fgetl(fid); end
 while ~feof(fid)
  str=fgetl(fid);
  [ tmp dg_t ]	= strread(str, '%f %f');
  dg		= [ dg	; dg_t ];
 end
 
 fclose(fid);				% close the dynamic file

 point			= DLS.Point;
 point.Instrument	= Instruments.ALV;
 point.T		= T;
 point.Angle		= angle;
 point.Tau_raw		= tau;
 point.G_raw		= g;
 point.dG_raw		= dg;		% this triggers the event in Point
 point.correct_G;

end

