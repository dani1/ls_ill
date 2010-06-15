function instrument = find_instrument( path );
% This function tries to understand what kind of instrument you have used according
% to some properties of your path

 % possible instruments
 instrument1	= 'Malvern Zetasizer Nano';
 instrument2	= 'ALV CGS3 and 7004/FAST';

 % PATHS: these can be updated with new versions of the program
 SLS_PATH_MALVERN	= [ path, '_SLS.txt' ];
 DLS_PATH_MALVERN	= [ path, '_DLS.txt' ];
 SLS_PATH_ALV		= [ path, '_bak.txt' ];
 for i = 1 : 100
  DLS_PATH_ALV{i}	= [ path, num2str(i-1,'%4.4u'), '.ASC'];
 end

 assignin('base','aa',SLS_PATH_ALV);

 % check for the existence of the right files
 if exist(SLS_PATH_MALVERN) == 2 || exist(DLS_PATH_MALVERN) == 2
  instrument	= instrument1;

 else
  instrument	= instrument2;

 end

end	% find_instrument
