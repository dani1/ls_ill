%============================================================================
%
% Utils class for common utility functions for the LS classes
%
% Author:	Fabio Zanini
% Companies:	Institut Laue Langevin, Universität Tübingen
% Date:		15 May 2010
% Version:	1.3
% 
% The idea is to collect some functions useful across all LS classes,
% which become subclasses of this one and inherit the methods
%
%============================================================================

classdef Utils < handle
 %===========================================================================
 % METHODS
 %===========================================================================
 methods ( Access = public )							% subclasses inherit the methods

  %==========================================================================
  % MONITOR PROPERTIES
  %==========================================================================
  function monitorprop ( self, varargin )
  % monitor selected properties through listeners

   props	= varargin;

   % add listeners (use anonymous functions...!)
   for i = 1 : length(props)
    try addlistener(self,props{i},'PostSet',@(src,evnt)prop_list(self,src,evnt) ); end	% monitor
   end

  end	% monitorprop

  %==========================================================================
  % ERASE DATAPOINT
  %==========================================================================
  function erase_datapoint (self, varargin )
  % erase one datapoint from Data , Gamma, etc. by specifying some arguments like
  % concentration, Q, etc.

   try	assert(nargin > 1)
   catch	error('Yo do not want to erase all data... choose C, Q and/or Angles.');
   end
 
   try	args	= struct(varargin{:});					% get arguments
   catch error('Arguments not understood!');
   end
 
   argnames	= fieldnames(args);
   len_before	= length(self.Data.C);					% the length of Data fields before cutting

   index	= logical(zeros(1,len_before));				% initialize that no datapoint survives
   for i = 1 : length(argnames)						% for every chosen argument...
    index = index | (self.Data.(argnames{i}) ~= args.(argnames{i}));	% ...let the elements of Data with a different value of that argument survive
   end

   self.Data	= self.correct_Data_fields ( index );			% correct the fields in Data
   self.correct_Datalike_properties(index, len_before);			% correct Gammas, etc.

  end	% erase_datapoint

  %==========================================================================
  % PROPERTIES LISTENERS
  %==========================================================================
  function prop_list (self, src, evnt )
  % if one of the observed fields has changed, update the Data struct and the other observed
  % properties to the new situation

   name		= src.Name;					% name of the property we monitor
   len_before	= length(self.Data.(name));			% the length of Data fields before cutting

   index	= ismember(self.Data.(name) , self.(name) );	% set to 1 if the respective property is found

   self.Data	= self.correct_Data_fields ( index );		% correct the fields in Data
   self.correct_Datalike_properties(index, len_before);		% correct Gammas, etc.

   self.update_uniques('C','Phi','Q','Q2','Angles');		% update unique props

  end	% prop_list

  %==========================================================================
  % CORRECT DATA FIELDS
  %==========================================================================
  function Data = correct_Data_fields (self, index )
  % correct the Data struct after events
  
   fnames	= fieldnames(self.Data);			% get all Data fields
   Data		= self.Data;					% copy the vector (for clarity)

   for i = 1 : length(fnames)					% for every field in Data...	
    Data.(fnames{i}) = self.Data.(fnames{i}) ( index );		% ...keep only the good data
   end

  end	% correct_data_fields


  %==========================================================================
  % CORRECT DATA-LIKE PROPERTIES
  %==========================================================================
  function correct_Datalike_properties (self, index, len_before )	% mmm, I would prefer a more transparent variant
  % correct the Data-like properties, such as the Gammas, after events

   props	= properties(self);				% properties of the calling class

   for i = 1 : length(props)					% scan all properties and...
    try
     assert( length(self.(props{i})) == len_before );		% ...check if the property is as long as a Data field...
     self.(props{i})	= self.(props{i})(index);		% ...if it is, reduce as for the Data struct.
    end
   end

  end	% correct_Datalike_properties

  %==========================================================================
  % UPDATE UNIQUE PROPERTIES
  %==========================================================================
  function update_uniques (self, varargin )			% mmm, I would prefer a more transparent variant
  % update unique properties such as C, Q, Angles

   argnames		= varargin;				% get the arguments

   for i = 1 : length(argnames)
    self.(argnames{i})	= unique(self.Data.(argnames{i}));	% recalculate uniques
   end

  end	% update_uniques

  %==========================================================================
  % UPDATE OR ADD DYNAMIC PROPERTIES
  %==========================================================================
  function check_add_prop ( self, varargin )
  % check whether a certain property is already in the class or not,
  % and add it if requested

   try		props = struct(varargin{:});			% get the properties
   catch
    error('There is some problem with your input syntax.');
   end
 
   f = fieldnames(props);					% find out the fieldnames

   for i = 1 : length(f)
    if ~ismember(f{i},properties(self))				% add the prop if necessary
     self.addprop(f{i});
    end
    self.(f{i}) = props.(f{i});					% fill the prop
   end

  end	% check_add_prop

  %==========================================================================
  % SAVE VARIABLES AS A COLUMN-SEPARATED ASCII FILE
  %==========================================================================
  function save_ascii ( self, filename, varargin )
  % save some variables as a matrix (column-separated) in an ascii text file
  % N.B.: the input vars MUST be column vectors!
  %
  % The file begins with a header line, describing the vars

   args		= struct(varargin{:});				% get the variables
   M		= [];						% initialize the matrix
   argnames	= fieldnames(args);				% initialize the name array
   filepath	= ['~/tm/LS/results/',filename];		% default folder TODO: choose a sensible one!
   headerline	= '#';


   for i = 1 : length(argnames)
    args.(argnames{i}) = args.(argnames{i})(:);			% they must be columns
    M		= [ M		args.(argnames{i})	];	% fill the matrix
    headerline	= [ headerline	'\t' argnames{i} 	];	% fill the header line
   end

   headerline = [headerline, '\n'];				% go to newline

   fid	= fopen(filepath,'w');					% open file
    fprintf(fid,headerline);					% Header line
   fclose(fid);

   save(filepath,'M','-ascii','-append');			% save the data

  end	% save ascii

  %============================================================================
  % SET PARAMETER FROM INDEPENDENT
  %============================================================================
  function parameter = set_parameter ( self, independent )

   switch independent
    case 'C'							% concentration
     parameter = 'Q2';

    case 'Phi'							% concentration
     parameter = 'Q2';

    case 'I'							% ionic strength
     parameter = 'Q2';

    case 'Q2'							% Q squared
     parameter = 'C';

    otherwise							% not recognized?!
     error('Independent not recognized!');
   end

  end	% set_parameter

 end	% methods

 %===========================================================================
 % STATIC METHODS
 %===========================================================================
 methods (Static)

  %============================================================================
  % AVERAGE
  %============================================================================
  function [ xm ym dym ] = average ( x, y, dy )
  % this function calculates the weighed average of a vector y using a
  % projection onto unique(x)

   if nargin == 2	dy = ones(1,length(y));		end				% non-weighed mean

   xm	= unique(x);									% decide the projection
   weights	= 1 ./ dy.^2;								% calculate weigths for the weighed sum
   for i = 1 : length(xm)
    index	= ( x == xm(i) );							% get the index...
    ym(i)	= sum ( y(index) .* weights(index) )	/ sum ( weights(index) );	% ...calculate y averaged...
    dym(i)	= sum ( dy(index) .* weights(index) )	/ sum ( weights(index) );	% ...calculate dy averaged
   end

  end	% average
  
 end	% static methods

end	% Utils class
