classdef Data < hgsetget	% necessary to inherit properties from the methods
%============================================================================
%
% Data class for loading SLS and DLS data
%
% Author:	Fabio Zanini
% Companies:	Institut Laue Langevin, Universität Tübingen
% Date:		20 Mar 2010
% Version:	1.1
%
% Syntax: dataclass = Data ( path, varargin )
%
% - path is the files location, like '/home/.../BSA_200_15_sample'. A relative
%   address can be used, but please mind your directory tree;
%
%   N.B.: if your instrument is the Malvern Zetasizer Nano, your filenames must be
%         the path + '_SLS.txt' for the static file and path + '_DLS.txt' for the
%         dynamic one.
%         If your instrument is the ALV, your filenames must be path + 0000.ASC
%         or any other startnumber and following for the dynamic files, and
%         path + '_bak.txt' for the static file!
%
% - varargin is a syntax of MatLab which allows for a variable number of input
%   options. The syntax is the following:
%
%        LSData ( path, 'option1', value1, 'option2', value2, ... )
%
%   Allowed options are at this moment the following:
%   - C :		real concentration (from UV; SLS/DLS)
%   - C_set :		concentration you set in the LS instrument (SLS)
%   - n :		real index of refraction of the solvent (from literature or
%			recfractometry;	SLS/DLS, is used to compute Q!)
%   ( n_set (NO!)	NOT AVAILABLE, read from raw data file			)
%   - dndc :		increment of index of refraction with sample concentration
%   - dndc_set :	increment of index of refraction set in the LS instrument
%   - Repetitions :	number of good repetitions at every angle (SLS/DLS)
%   - Runstart :	number of the first DLS file (e.g. DLSfile_0002.ASC). This
%			is passed directly to a low-level routine, so you can leave
%			it away if you do not need it for your instrument
%   - Runs :		the number of DLS data files. If you let it empty, a low-
%			level routine is called to find it out. However, since this
%			routine is instrument-specific, you have to write it if you
%			are not measuring with the ALV CGS3 at ILL!
%
% ---- Rationale of the class ----
%
% When you open an instance of Data, its constructor is called. It performs  the 
% following tasks:
%
% - eats the path and stores it
% - tries to understand the optional arguments, and sets the default values for
%   those let empty by the user
% - loads general information about the data, e.g. instrument name, wavelength of
%   the laser, etc. This is done through a low-level function. If you plan to use
%   this class outside the ILL, you'll need to write it again
% - loads static data. This is done at high-level inside this class, and at low-
%   level through an external, instrument-specific function. See the note above
% - corrects static data by asking the user or falling back to some default
%   corrections. It is cumbersome to bother the user every time, but it is really
%   important for him to realize that they are needed and sample-specific
% - calculates Q vectors for the static data (using the aforementioned corrections)
% - loads the dynamic data. Again, this is done by a joint internal high-level and
%   external low-level pair of functions. Since DLS is a very sensitive technique,
%   some kind of data filtering is necessary at this point (unless you want to
%   see beautiful dust in most of your datasets). This may be performed
%   automatically by the high-level function, but the filtering criteria are
%   somehow arbitrary. Otherwise, a manual selection can be activated. Please,
%   take a look at that function.
% - calculates Q vectors for the dynamic data (these need not be the same as before)
%
% ----- SLS -----
%
% It is necessary to correct the SLS data because of some parameters like
% the sample concentration, which are unknown at measurement time. They are
% the following:
%
% - refractive index n (used only to calculate q)
% - refractive index increment dn/dc
% - concentration as measured with UV absorption
%
% Starting from version 1.1, the class tries to be less user-invasive. If the user
% states some of these three parameters as optional args, they are taken in. Otherwise,
% the class tries a standard data correction according to "normal conditions", and
% tells the user about this.
%
% At present, it is not possible to modify the parameters at runtime. This creates
% quite a few problems, like double-corrections, user bothering, and the need of
% event-aware properties. Maybe in the future...
%
% ----- DLS -----
%
% DLS error bars on the Gs are calculated *not* as stddev between the
% measurement repetitions, beacuse one would need tens of repetitions.
% Rather, stddev of G is calculated according to Schaetzel and Peters as of:
% W. Brown (ed.), Dynamic Light Scattering, Clarendon Press, Oxford, 1993.
% This should not be a problem since later-induced statistical incertainties
% tend to dominate the final error on the diffusion constant.
%
%============================================================================

 %============================================================================
 % PUBLIC PROPERTIES
 %============================================================================
 properties ( Access = public )
  Path
  % path is the path of the correlation files, w/o 0XXX and w/ extension
  % Ex: /home/.../BSA_200_15_sample

  % general
  Instrument
  Lambda
  Unit_Lambda
  Unit_Q	= 'A^{-1}';	% Unit_Q is always A^-1
  C				% read from varargin or asked
  C_set				% read from varargin or asked
  Unit_C	= 'mg/ml';
  n				% asked or n_set
  n_set
  dndc		= 0.1850;	% usual value in Tuebingen
  dndc_set	= 0.1820;	% usual value at ILL
  Unit_dndc	= 'ml/g'

  % SLS
  Angles_static
  Q_static
  KcR
  dKcR
  Unit_KcR	= 'mol/g';	% standard units

  % DLS
  Repetitions
  Angles_dyn
  Q_dyn
  Tau
  Unit_Tau	= 'ms';		% value at ILL
  G
  dG

 end	% public_properties

 %============================================================================
 % PRIVATE PROPERTIES
 %============================================================================
 properties ( Access = private )

  Runstart	= 0;
  Runs

 end	% private properties

 %============================================================================
 % PUBLIC METHODS
 %============================================================================
 methods ( Access=public )

  %============================================================================
  % CONSTRUCTOR
  %============================================================================
  function lsd = Data ( path, varargin )

   lsd.Path = path;			% path
   lsd.set_optional ( varargin );	% optional arguments like c_set, ecc.

   % the first thing is to get the general information about the instrument, lambda, etc.
   [ lsd.Instrument lsd.Lambda lsd.Unit_Lambda lsd.n_set ] = lsd.load_general_data ( varargin );

   [ lsd.Angles_static lsd.KcR lsd.dKcR ] = lsd.load_static_data ( varargin );	% SLS data (corrected)
   lsd.Q_static = lsd.compute_Q_from_angles ( lsd.Angles_static );		% calculate Q

   lsd.load_dyn_data;								% DLS data and correction (auto or manual)
   lsd.Q_dyn = lsd.compute_Q_from_angles ( lsd.Angles_dyn );			% calculate Q

  end	% constructor 

 end	% public methods

 %============================================================================
 % PRIVATE METHODS
 %============================================================================
 methods ( Access = private )

  %=========================================================================================
  % SET OPTIONAL PARAMETERS OF THEIR DEFAULT VALUES
  %=========================================================================================
  function set_optional ( obj, argvect )
  % set optional class properties read from varargin of the class itself

   if mod(length(argvect),2)
    error('You should choose parameter names and input a values for them.');
   else
    options		= struct(argvect{:});				% read the optional args

    % Repetitions
    try
     obj.Repetitions	= options.Repetitions;
    end

    % Runstart
    try
     obj.Runstart	= options.Runstart;
    end

    % Instrument
    try
     obj.Instrument	= options.Instrument;
    end

   end

  end	% set_optional

  %=========================================================================================
  % LOAD GENERAL DATA
  %=========================================================================================
  function [ instrument lambda unit_lambda n_set ] = load_general_data ( obj, argvect )
  % this is just a switch function. We have to understand which kind of instrument
  % we're using, and then call external functions accordingly

   PATH_ALV		= [ obj.Path, '0', num2str(obj.Runstart,'%3.3u'), '.ASC'];			% this is the syntax for ALV files
   PATH_MALVERN_DLS	= [ obj.Path, '_DLS.txt'];							% this is the syntax for Malvern
   PATH_MALVERN_SLS	= [ obj.Path, '_SLS.txt'];							% this is the syntax for Malvern
 
   if ( exist(PATH_ALV) == 2	|| ( ~isempty(obj.Instrument) && regexp(obj.Instrument,'ALV') ) )				% ALV?
    instrument	= 'ALV CGS3 and 7004/FAST';
    [ lambda unit_lambda n_set ] = obj.read_general_file_ALV ( obj.Path, obj.Runstart );

    % Runs
    try
     obj.Runs	= options.Runs;
    catch
     obj.Runs	= obj.find_runs(obj.Path,obj.Runstart);			% find the Runs if necessary (may be incorrect)
    end
 
   elseif ( exist(PATH_MALVERN_SLS) == 2 || exist(PATH_MALVERN_DLS) == 2 || ( ~isempty(obj.Instrument) && regexp(obj.Instrument,'Malvern') ) )			% Malvern?
    instrument	= 'Malvern Zetasizer Nano';
    [ lambda unit_lambda n_set ]	= obj.set_general_parameters_Malvern;

    % Runs
    try
     obj.Runs	= options.Runs;
    catch
     obj.Runs	= 1;			% the standard for the Malvern instrument is 1
    end

 
   else												% none of them?!
    fprintf('\nEh? Instrument not recognized. Are you sure that your path is correct?\n\n');
   end

  end	% load_general_data

  %=========================================================================================
  % LOAD STATIC DATA
  %=========================================================================================
  function [ angles_static kcr dkcr ] = load_static_data ( obj, argvect )
  % reading static data from text files and correcting for them at initialization
  
   if regexp(obj.Instrument,'Malvern')					% Malvern instrument?

    PATH_MALVERN_SLS	= [ obj.Path, '_SLS.txt'];					% this is the syntax for Malvern

    [ angles_static kcr dkcr ] = obj.read_static_file_Malvern( PATH_MALVERN_SLS );

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
 

   elseif regexp(obj.Instrument,'ALV')					% ALV Instrument?
    [ Ang_a KcR_a dKcR_a ] = obj.read_static_file_ALV ( obj.Path );	% read raw data (with repetitions)

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
    options	= struct(argvect{:});					% read the optional args
  
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

   else
    fprintf('\nInstrument not recognized!\n\n');
   end
   
  end	% load_static_data 

  %=========================================================================================
  % READ CORRELATION DATA
  %=========================================================================================
  function load_dyn_data ( obj )
  % This one is tricky. You have to find all the correlation files, find out the number
  % of repetitions. Every time that for one angle there are more gs, plot them and ask
  % the user what to do. If the user does not say anything, choose the gs with the smallest
  % integral in the intermediate lagtime range (fewest aggregates).

   if regexp(obj.Instrument,'Malvern')					% Malvern instrument?

   PATH_MALVERN_DLS	= [ obj.Path, '_DLS.txt'];			% this is the syntax for Malvern

    [ Ang_a Tau_a G_a dG_a ] = obj.read_dyn_files_Malvern( PATH_MALVERN_DLS );
    obj.Repetitions	= 1;						% on Malvern, we make a mean of all measurements


   elseif regexp(obj.Instrument,'ALV')					% ALV Instrument?

    [ Ang_a Tau_a G_a dG_a ] = obj.read_dyn_files_ALV( obj.Path, obj.Runstart, obj.Runs );	% read raw data
  
   end
  
   Ang_unique = unique(Ang_a);				% unique vector of angles

   % find out the number of repetitions IF the parameter was not passed explicitely!
   if isempty(obj.Repetitions)
    for i = 1 : obj.Runs
     rep(i) = sum(ismember(Ang_a,Ang_a(i)));
    end
    obj.Repetitions	= min(rep);
   end
  
   % fill the Angles vector using the Repetitions
   for i = 1 : length(Ang_unique)
    for k = 1 : obj.Repetitions
     obj.Angles_dyn(obj.Repetitions*(i-1)+k)	= Ang_unique(i);
    end
   end
  
   % pick out the best data points asking the user. the default choice is based on
   % some quality parameter ( smallest integral in a definite range, etc. ).
   for i = 1 : length(Ang_unique)
  
    % create an empty struct with the candidate vectors
    candidates	= struct('tau',[],'g',[],'dg',[]);
    for j = 1 : length(Ang_a)
     if Ang_a(j) == Ang_unique(i)
      candidates.tau	= [ candidates.tau	{Tau_a{j}} ];
      candidates.g	= [ candidates.g	{G_a{j}} ];
      candidates.dg	= [ candidates.dg	{dG_a{j}} ];
     end
    end

    % normalize candidates
    for j = 1 : length(candidates.tau)
     normali(j)		= mean_correlation( candidates.g{j}, 1, 20 );
     candidates.g{j}	= candidates.g{j} ./ normali(j);
     candidates.dg{j}	= candidates.dg{j} ./ normali(j);
     clear normali;
    end
  
    % if the number of candidates is equal to the number of places, do not show anything
    if length(candidates.tau) == obj.Repetitions

     for k = 1 : obj.Repetitions
      obj.Tau{obj.Repetitions*(i-1)+k}	= candidates.tau{k};
      obj.G{obj.Repetitions*(i-1)+k}	= candidates.g{k};
      obj.dG{obj.Repetitions*(i-1)+k}	= candidates.dg{k};
     end
  
    % otherwise, plot the candidates with their (new) index and filter them
    else
%     figure;									% open figure
%     hold all;
%     set( gcf, 'color', 'white' );
%     set( gca, 'box', 'on', 'xscale','log' );
%     xlabel(['\tau [ ms ]']);
%     ylabel('g_2-1/\beta [ normalised ]');
%     pl		= [];
%     leg		= [];
%   
%     for j = 1 : length(candidates.tau)
%      pl		= [ pl;  plot(candidates.tau{j},candidates.g{j}) ];	% plot
%      leg		= [ leg; ['Plot n.',num2str(j,'%2.2u')] ];		% legend
%     end
%     legend(pl,leg);
  
     % calculate the best functions based on:
%     index = find_minimal_integrals ( candidates, obj.Repetitions );		% the minimal integral criterium
%     index = find_minimal_set ( candidates, obj.Repetitions );			% the minimal value at some points
     index = find_fastest_thre ( candidates, obj.Repetitions );			% the fastest decay to some g threshold
  
     % ask the user to choose some correlation function
%     answer = input(['Please choose ', num2str(obj.Repetitions),' correlations: [ ',num2str(index),' ] ']);
%     if isempty(answer)
      answer = index;
%     end
  
%     close(gcf);								% close figure
  
     % insert the chosen correlations
     for k = 1 : obj.Repetitions
      obj.Tau{obj.Repetitions*(i-1)+k}	= candidates.tau{answer(k)};
      obj.G{obj.Repetitions*(i-1)+k}	= candidates.g{answer(k)};
      obj.dG{obj.Repetitions*(i-1)+k}	= candidates.dg{answer(k)};
     end
 
    end
  
   end


  end	% load_dyn_data 

  %=========================================================================================
  % COMPUTE Q FROM ANGLES
  %=========================================================================================
  function q = compute_Q_from_angles ( obj, angles )
  % fill the Q vector with lsd.Runs entries starting from the angles (both SLS and DLS)

   % check the unit conversion ( convert lambda into angstrom )
   if strcmpi(obj.Unit_Lambda,'A')		lambda = obj.Lambda;
   elseif strcmpi(obj.Unit_Lambda,'nm')		lambda = 10 * obj.Lambda;
   else				warning('Wavelength unit not supported.');
   end

    q = 4 * pi * obj.n * sind( 0.5 * angles ) ./ lambda;	% compute the q vectors
  end	% compute_Q_from_angles 

 end	% private methods

 %============================================================================
 % STATIC METHODS (private)
 %============================================================================
 methods ( Access = private, Static )

  %=========================================================================================
  % find number of runs (separate file)
  %=========================================================================================
  runs = find_runs (path, runstart);

  %=========================================================================================
  % read_general_file (separate file)
  %=========================================================================================
  [ lambda unit_lambda n_set ] = read_general_file_ALV ( path, runstart );
  [ lambda unit_lambda n_set ] = set_general_parameters_Malvern;

  %=========================================================================================
  % read_static_file (separate file)
  %=========================================================================================
  [ ang kcr dkcr ] = read_static_file_ALV ( path );
  [ ang kcr dkcr ] = read_static_file_Malvern ( path );

  %=========================================================================================
  % read_dyn_files (separate file)
  %=========================================================================================
  [ angle tau g dg ] = read_dyn_files_ALV ( path, runstart, runs );
  [ angle tau g dg ] = read_dyn_files_Malvern ( path );

 end	% static methods

% end of class
end

%=========================================================================================
%
% LOW LEVEL ROUTINES
%
%=========================================================================================
%=========================================================================================
% filter the static data according to minimal relative standard deviation
%=========================================================================================
function [ kcr dkcr ] = filter_stddev ( candidates )
 [ dkcr index ]	= min(candidates.dkcr);
 kcr		= candidates.kcr(index);
end	% filter_std_dev

%=========================================================================================
% find mean correlation function between two thresholds
%=========================================================================================
function i_c	= mean_correlation( g, threshold1, threshold2)
 i_c	= 0;
 for i = threshold1 : threshold2
  i_c = i_c + g(i) / ( threshold2 - threshold1 );
 end
end	% mean_correlation

%=========================================================================================
% calculate index for minimum integral of correlations between two thresholds
%=========================================================================================
function index = find_minimal_integrals ( cand, rep )

 thresh_min	= 0.8;	thresh_max	= 0.2;	% set the threshold for g

 % find the indices of the thresholds
 ind_min = ones(length(cand.tau),1);		ind_max = ones(length(cand.tau),1);
 for i = 1 : length(cand.tau)
  tot = [ 1: length(cand.tau{i}) ];
   ind_min(i) = max(tot(cand.g{i} > thresh_min));
   ind_max(i) = max(tot(cand.g{i} > thresh_max));
 end

 % calculate the mean values (integrals)
 for i = 1 : length(cand.tau)
  mean_corr(i)	= mean_correlation ( cand.g{i}, ind_min(i), ind_max(i) );
 end
 for k = 1 : rep
  [ m index(k) ]	= min(mean_corr);	% write the index with the minimum integral
  mean_corr(index(k))	= max(mean_corr)+1;	% exclude that value from the next minumum
 end
end	% find_minimal_integrals

%=========================================================================================
% calculate index for minimum correlation at some points
%=========================================================================================
function index = find_minimal_set ( cand, rep )

 set = [ 1e-2 1e-1 1 1e1 1e2 ];

 % find the indices of the set points and calculate g at those points
 ind = ones(length(cand.tau),length(set));
 for i = 1 : length(cand.tau)
  tot = [ 1: length(cand.tau{i}) ];
  for j = 1 : length(set)
   ind(i,j) = max(tot(cand.tau{i} < set(j)));
   values(i,j)	= cand.g{i}(ind(i,j));
  end
  weighed(i) = sum(values(i,:));		% weight the points
 end

 for k = 1 : rep
  [ m index(k) ]	= min(weighed);		% write the index with the minimum weight
  weighed(index(k))	= max(weighed)+1;	% exclude that value from the next minumum
 end

end	% find_minimal_set

%=========================================================================================
% calculate index for fastest decay to some g threshold
%=========================================================================================
function index = find_fastest_thre ( cand, rep )

 thre = 2e-2;

 for i = 1 : length(cand.tau)
  value(i) = max(cand.tau{i}(cand.g{i} > thre )); % find the maximal time above the thre
 end

 for k = 1 : rep
  [ m index(k) ]	= min(value);		% write the index with the minimum weight
  value(index(k))	= max(value)+1;		% exclude that value from the next minumum
 end

end	% find_fastest_thre

%=========================================================================================
% calculate minimum of time steps of an array of correlations
%=========================================================================================
function last_time = min_last_time ( tau )
 for i = 1 : length(tau)
  l(i)	= length(tau{i});
 end
 last_time	= min(l);
end	% min_last_time
