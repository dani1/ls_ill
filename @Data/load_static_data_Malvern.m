  function [ angles_static kcr dkcr ] = load_static_data_Malvern ( obj, slspath, options )

   [ angles_static kcr dkcr ] = obj.read_static_file_Malvern( slspath );	% function who read the files low-level

   % Now come some corrections

   % C (no fallback)
   try
    obj.C	= options.C;
    catch
    obj.C	= input('What is the sample concentration after UV [mg/ml]? [ no default ] ');
    if isempty(obj.C)
     error('You HAVE to insert some concentration value, how are you to interpret the data without it?');
    end
   end

   % n (fallback: n_set, from read_general_file_* )
   try
    obj.n	= options.n;
   catch
    obj.n	= obj.n_set;
    fprintf('No corrected index of refraction for the solvent has been input. Falling back to n_set\n');
   end
 
  end	% load_static_data_Malvern
