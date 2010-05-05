function [ ang kcr dkcr ] = read_static_file_ALV ( path )
% read SLS data from files for the ALV CGS3 goniometer at ILL, Grenoble

  % initialize to void, if the function outputs void then skip SLS
  ang	= [];
  kcr	= [];
  dkcr	= [];
 
 static_file= [path,'_bak.txt'];				% create static file
 if exist(static_file,'file') ~= 2				% check for existance
  static_file = input([static_file, ...
	' not found. Please enter the full path of the static file or let empty for skipping: '],...
			's');
 end

 % check whether the new file exists
 if exist(static_file,'file') == 2

  fid = fopen(static_file);	% open the file
 
  while ~feof(fid)
   str = fgetl(fid);
   if str(1) ~= 'A'		% check whether this is a header or a useful line

    % get the data about angles, kcr and dkcr
    [ ang_t tmp tmp tmp tmp tmp tmp tmp kcr_t dkcr_t tmp tmp ] = strread(str, '%f %f %f %f %f %f %f %f %f %f %s %s');

    ang				= [ ang ang_t ];
    kcr				= [ kcr kcr_t ];
    dkcr			= [ dkcr dkcr_t ];

   end
  end
 
  fclose(fid);			% close the file

 end

end	% read_static_file
