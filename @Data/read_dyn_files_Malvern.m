function [ Angles Tau G dG ] =  read_dyn_files_Malvern ( path )
% read text files to get DLS data
% aaaaah, the Malvern output format is mad mad mad mad!
 dynamic_file	= path;


 t	= [];
 g	= [];
 dg	= [];

 if exist(dynamic_file) ~= 2					% check file exsistence

  error('File not found!');

 else

  fid = fopen(dynamic_file);					% open the file

  while ~feof(fid)						% till the end of file, do the following

   str = fgetl(fid);						% get the lines

   if ( ~isempty(regexp(str,'Correlation Delay Times')) & ~isempty(regexp(str,'Correlation Data')) )	% get the line with the format

    formatvect	= textscan(str,'%s','delimiter', '\t');		% convert into a cell array
    formatvect	= formatvect{1};

    str = fgetl(fid);						% this line has the data

    datavect	= textscan(str,'%s','delimiter', '\t');		% convert data sting to cell array (there is also other info)
    datavect	= datavect{1};


    for i = 1 : length(datavect)

     if		~isempty(regexp(formatvect{i},'Correlation Delay Times'))	% add lag times
      t = [ t str2num(regexprep(datavect{i},',','.')) ];
     elseif	~isempty(regexp(formatvect{i},'Correlation Data'))		% add correlogramm
      g = [ g str2num(regexprep(datavect{i},',','.')) ];
     end
    end

   end


  end

  fclose(fid);

 end


 % store the info
 Angles	= 173;
 Tau	= {t};
 G	= {g};
 dG	= {0};

end	% read_dyn_files_Malvern
