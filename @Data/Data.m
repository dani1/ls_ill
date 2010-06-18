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
%   - T :		the temperature of the measurement
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
  GENPath
  % GENPath is the path where the general properties are stored

  SLSPath
  % SLSpath is the path of the table, inclusive the extension

  DLSPath
  % DLSpath is the path of the correlation files, w/o 0XXX and w/ extension
  % Ex: /home/.../BSA_200_15_sample

  % general
  Sample	= 'BSA'		% default sample: Bovine Serum Albumin
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
  T		= 22-LIT.Constants.T0;	% temperature (default: 22°C)

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
  function lsd = Data ( varargin )

   try		args	= struct(varargin{:});
   catch	error('Please check your syntax.');	end			% try to read the arguments

   try		lsd.Sample	= args.Sample;
   catch	disp('Default sample: BSA.');		end			% Sample

   try		lsd.Repetitions	= args.Repetitions;	end			% Repetitions
   try		lsd.Runstart	= args.Runstart;	end			% Runstart

   try		lsd.C	= args.C;	end			% Runstart
   try		lsd.C_set	= args.C_set;	end			% Runstart
   try		lsd.n	= args.n;	end			% Runstart

   try		lsd.T	= args.T;						% Temperature
    		if	lsd.T < 150	lsd.T = lsd.T-LIT.Constants.T0;	end	% check for °C input temperatures
   catch	disp(	['T =' num2str(lsd.T+LIT.Constants.T0),'°C.']	);	% fall back to 22°C
   end

   try		lsd.DLSPath	= args.DLS;
   catch	disp('No DLS data.');			end			% DLS

   try		lsd.GENPath	= args.GEN;					% General parameters
   catch
    try		lsd.GENPath	= args.DLS;
    catch	error('No path for the general parameters!');
    end
   end

   try		[ lsd.Instrument f ]	= lsd.find_instr( args.Instrument );	% Instrument
   catch	[ lsd.Instrument f ]	= lsd.find_instr;	end

   [	lsd.Lambda, 		...
	lsd.Unit_Lambda,	...
	lsd.n_set 	]	= lsd.(f.gen) ( lsd.GENPath,lsd.Runstart );	% general options

   try										% SLS
    lsd.SLSPath	= args.SLS;
    [	lsd.Angles_static,	...
	lsd.KcR,		...
	lsd.dKcR 	]	= lsd.(f.sta) ( lsd.SLSPath, args 	);	% load SLS data
    lsd.Q_static		= lsd.Q_from_angles ( lsd.Angles_static );	% calculate Q

   catch	disp('No SLS data.');			end


   try		lsd.Runs = options.Runs;					% Runs
   catch	lsd.Runs = lsd.(f.fin) ( lsd.DLSPath, lsd.Runstart );   end

   try										% DLS if needed
    assert(~isempty(lsd.DLSPath));
    [	Ang_a,	...
	Tau_a,	...
	G_a,	...
	dG_a 	] = lsd.(f.dyn) ( lsd.DLSPath, lsd.Runstart, lsd.Runs 	);	% load DLS data
    [	lsd.Angles_dyn,	...
	lsd.Tau,	...
	lsd.G,		...
	lsd.dG 		] = lsd.filter_dyn_data(Ang_a,Tau_a,G_a,dG_a	);	% filter DLS data
    lsd.Q_dyn	= lsd.Q_from_angles ( lsd.Angles_dyn );			% calculate Q
   end

  end	% constructor 

 end	% public methods

 %============================================================================
 % PRIVATE METHODS
 %============================================================================
 methods ( Access = private )

  %=========================================================================================
  % FILTER CORRELATION DATA
  %=========================================================================================
  function [ Angles Tau G dG ] = filter_dyn_data ( obj, Ang_a, Tau_a, G_a, dG_a, options )
  % This one is tricky. You have to find all the correlation files, find out the number
  % of repetitions. Every time that for one angle there are more gs, plot them and ask
  % the user what to do. If the user does not say anything, choose the gs with the smallest
  % integral in the intermediate lagtime range (fewest aggregates).

   % find repetitions
   if isempty(obj.Repetitions);
    obj.Repetitions	= find_repetitions(Ang_a); 
   end

   Ang_unique = unique(Ang_a);				% unique vector of angles
 
   % fill the Angles vector using the Repetitions
   for i = 1 : length(Ang_unique)
    for k = 1 : obj.Repetitions
     Angles(obj.Repetitions*(i-1)+k)	= Ang_unique(i);
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
      Tau{obj.Repetitions*(i-1)+k}	= candidates.tau{k};
      G{obj.Repetitions*(i-1)+k}	= candidates.g{k};
      dG{obj.Repetitions*(i-1)+k}	= candidates.dg{k};
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
      Tau{obj.Repetitions*(i-1)+k}	= candidates.tau{answer(k)};
      G{obj.Repetitions*(i-1)+k}	= candidates.g{answer(k)};
      dG{obj.Repetitions*(i-1)+k}	= candidates.dg{answer(k)};
     end
 
    end
  
   end


  end	% load_dyn_data 

  %=========================================================================================
  % COMPUTE Q FROM ANGLES
  %=========================================================================================
  function q = Q_from_angles ( obj, angles )
  % fill the Q vector with lsd.Runs entries starting from the angles (both SLS and DLS)

   % check the unit conversion ( convert lambda into angstrom )
   if strcmpi(obj.Unit_Lambda,'A')		lambda = obj.Lambda;
   elseif strcmpi(obj.Unit_Lambda,'nm')		lambda = 10 * obj.Lambda;
   else				warning('Wavelength unit not supported.');
   end

    q = 4 * pi * obj.n * sind( 0.5 * angles ) ./ lambda;	% compute the q vectors
  end	% compute_Q_from_angles 

  %=========================================================================================
  % LOAD STATIC DATA
  %=========================================================================================
  [ angles_static kcr dkcr ] = load_static_data_Malvern ( obj, slspath, options )
  [ angles_static kcr dkcr ] = load_static_data_ALV ( obj, slspath, options )

 end	% private methods

 %============================================================================
 % STATIC METHODS
 %============================================================================
 methods ( Access = private, Static )

  %=========================================================================================
  % find instrument
  %=========================================================================================
  function [ instrument functions ] = find_instr( instr );
  % This function tries to understand what kind of instrument you have used according
  % to some properties of your path

   % possible instruments
   Malvern	= 'Malvern Zetasizer Nano';
   ALV		= 'ALV CGS3 and 7004/FAST';

   if		nargin	== 0;		instrument = ALV;
   elseif	strfind(Malvern,instr)	instrument = Malvern;
   elseif	strfind(ALV,	instr)	instrument = ALV;
   else		error('Instrument not found.');
   end

   instr		= strtok(instrument);				% short version of the instrument

   % choose the methods according to the instrument
   functions.gen	= ['read_general_file_',instr];
   functions.sta	= ['load_static_data_',instr];
   functions.fin	= ['find_runs_',instr];
   functions.dyn	= ['read_dyn_files_',instr];
  
  end	% find_instrument 

  %=========================================================================================
  % find runs
  %=========================================================================================
  runs	= find_runs_ALV ( dlspath, runstart )
  runs	= find_runs_Malvern ( dlspath, runstart )

  %=========================================================================================
  % read_general_file
  %=========================================================================================
  [ lambda unit_lambda n_set ] = read_general_file_ALV ( genpath, runstart )
  [ lambda unit_lambda n_set ] = read_general_file_Malvern ( genpath, runstart )

  %=========================================================================================
  % read_static_file
  %=========================================================================================
  [ ang kcr dkcr ] = read_static_file_ALV ( slspath )
  [ ang kcr dkcr ] = read_static_file_Malvern ( slspath )

  %=========================================================================================
  % read_dyn_files
  %=========================================================================================
  [ angle tau g dg ] = read_dyn_files_ALV ( path, runstart, runs )
  [ angle tau g dg ] = read_dyn_files_Malvern ( path, runstart, runs )

 end	% static methods

% end of class
end

%=========================================================================================
%
% LOW LEVEL ROUTINES
%
%=========================================================================================

%=========================================================================================
% find the number of repetitions
%=========================================================================================
function repetitions = find_repetitions ( ang_a )
 for i = 1 : length(ang_a)
  rep(i) = sum( ang_a == ang_a(i) );		% calculate how often every element is repeated
 end
 repetitions	= min(rep);			% take the minimum
end	% find_repetitions
 
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
