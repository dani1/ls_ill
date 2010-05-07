function runs = find_runs_ALV ( dlspath, runstart );
% calculate the runs for a DLS measurement by the ALV instrument

 runs = 0;
 while exist([dlspath,num2str(runstart+runs,'%4.4u'),'.ASC'], 'file') == 2
  runs=runs+1;
 end

end	% find_runs_ALV
