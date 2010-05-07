%=========================================================================================
% If the user does not not pass the argument about the number of runs, try to find it out
% yourself, checking the number of files of the form:
%
% <sample_name>00XX.ASC
%=========================================================================================
function runs = find_runs (path, runstart, instrument)

 if ~isempty(regexp(instrument,'ALV'))

  runs = 0;
  while exist([path,'00',num2str(runstart+runs,'%2.2u'),'.ASC'], 'file') == 2
   runs=runs+1;
  end

 elseif ~isempty(regexp(instrument,'Malvern'))
  runs	= 1;
 end

end	% find_runs
