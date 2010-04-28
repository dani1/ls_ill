function [ files angles_one repetitions ] = load_files ( path, runstart, runs)

 % store the full path of all files into a temporary c.a.
 files_tmp	= cell (runs,1);
 angles		= zeros(runs,1);
 for j = 1 : runs
   jstr= num2str(j-1+runstart,'%04u');
   files_tmp{j} = [path,jstr,'.ASC'];
 end

 % read the files and store their angles
 for j = 1 : runs
  fid = fopen(files_tmp{j});
  while angles(j) == 0
   str = fgetl(fid);
   if strfind(str,'Angle')
    [ tmp tmp tmp angles(j) ] = strread(str, '%s %s %s %f');
   end
  end
  fclose(fid);
 end

 % count what is the minimal number of
 % repetitions and store it as every
 % (hopefully at least one angle is good!)
 angles_one 	= angles(1);
 counts		= 1;
 for j = 2 : runs
  if angles(j) ~= angles_one(length(angles_one))
   angles_one 	= [ angles_one angles(j) ];
   counts 	= [ counts 1 ];
  else
  counts(length(counts)) = counts(length(counts)) + 1;
  end
 end
 repetitions = min(counts);

 % fill the files vector with the right files
 % (the last ones)
 ss=0;
 files = cell(repetitions * length(angles_one) ,1 );
 for j = 1 : length(angles_one)
  ss=ss+counts(j);
  for sub_j = 1 : repetitions
   files{repetitions*(j-1)+sub_j} = files_tmp( ss-repetitions+sub_j );
  end
 end

 %end of function
end
