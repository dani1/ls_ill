function point = read_static_file (self, path )
% read SLS data from files for the ALV CGS3 goniometer at ILL, Grenoble
 i = 0;
 point     = SLS.Point;
 min_error = 1e-4;                              % smallest accepted relative error

 fid       = fopen( path );                           % open the file
 
 while ~feof(fid)

  str = fgetl(fid);

  try
   assert( str(1) ~= 'A' & ~isempty(str) );
   tmp   = textscan(str,'%f %f %f %f %f %f %f %f %f %f %s %s');
   angle = tmp{1};
   T     = tmp{7};
   kcr   = tmp{9};
   dkcr  = tmp{10};
   dkcr  = max(dkcr,min_error);

   i = i+1;
   point(i) = SLS.Point;
   point(i).Instrument = self;
   point(i).T          = T;
   point(i).Angle      = angle;
   point(i).KcR_raw    = kcr;
   point(i).dKcR_raw   = 0.01 * dkcr * kcr;                   % percentual errors!

  end

 end

 fclose(fid);           % close the file

end	% read_slspath
