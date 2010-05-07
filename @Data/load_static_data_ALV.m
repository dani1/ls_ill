function [ angles_static kcr dkcr ] = load_static_data_ALV ( obj, slspath, options )

 [ Ang_a KcR_a dKcR_a ] = obj.read_static_file_ALV ( slspath );	% read raw data (with repetitions)

 if isempty(Ang_a)
  fprintf('Skipping SLS data');					% print a message if skipping
 else

  angles_static	= unique(Ang_a);				% unique already orders the vector   

  % pick out the best data point on the basis of a quality parameter ( percentual error, smallest I, etc. );
  for i = 1 : length(angles_static)
   candidates = struct('kcr',[],'dkcr',[]);				% create empty candidate
   for j = 1 : length(Ang_a)
    if Ang_a(j) == angles_static(i)
     candidates.kcr	= [ candidates.kcr	KcR_a(j) ];
     candidates.dkcr	= [ candidates.dkcr	dKcR_a(j)];
    end
   end
  [ kcr(i) dkcr(i) ]	= filter_stddev ( candidates );			% filter stddev
  end   

  dkcr	= 0.01 .* dkcr .* kcr;					% turn relative into absolute errors
 end

 % Now we need to correct the parameters as specified by the user.
 % The protocol for every parameter is the following:
 % 1. look for the parameter in the optional args. If found, take it into the class
 % 2. If not found, use the fallback value of the class if present (warn the user)
 % 3. If no fallback value is found (e.g. for C_set), ask the user.

 % C_set (no fallback)
 try
  obj.C_set	= options.C_set;
 catch
  obj.C_set	= input('What is the concentration set in the LS instrument [mg/ml]? [ no default ] ');
  if isempty(obj.C_set)
   error('You HAVE to insert some concentration value, how are you to interpret the data without it?');
  end
 end

 % C (fallback: C_set)
 try
  obj.C	= options.C;
 catch
  obj.C	= obj.C_set;
  fprintf('No UV-determined concentration has been input. Falling back to C_set\n');
 end

 % n (fallback: n_set, from read_static_file_* )
 try
  obj.n	= options.n;
 catch
  obj.n	= obj.n_set;
  fprintf('No corrected index of refraction for the solvent has been input. Falling back to n_set\n');
 end

 % dndc (fallback: class default)
 try
  obj.dndc	= options.dndc;
 catch
  fprintf(['No corrected dn/dc has been input. Falling back to the standard: ', num2str(obj.dndc),' ml/g\n']);
 end

 % dndc_set (fallback: class default)
 try
  obj.dndc_set = options.dndc_set;
 catch
  fprintf(['No corrected dn/dc_set has been input. Falling back to the standard: ', num2str(obj.dndc_set),' ml/g\n\n']);
 end

 % correction of values
 correction_factor = ( ( obj.n * obj.dndc ) ^2 * obj.C ) ./ ( ( obj.n_set * obj.dndc_set ) ^2 * obj.C_set  );
 obj.KcR	= obj.KcR * correction_factor;
 obj.dKcR	= obj.dKcR* correction_factor;


 %=========================================================================================
 % NESTED FUNCTION: filter the static data according to minimal relative standard deviation
 %=========================================================================================
 function [ kcr dkcr ] = filter_stddev ( candidates )
  [ dkcr index ]	= min(candidates.dkcr);
  kcr		= candidates.kcr(index);
 end	% filter_std_dev

end	% load_static_data 
