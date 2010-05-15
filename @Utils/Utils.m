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
    try addlistener(self,props{i},'PostSet',@(src,evnt)prop_list(self,src,evnt) ); end	% monitor C
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

   self.update_uniques('Q','Angles');				% update unique props

  end	% prop_list

  %==========================================================================
  % CORRECT DATA FIELDS
  %==========================================================================
  function Data = correct_Data_fields (self, index )
  % correct the Data struct after events
  
   fnames	= fieldnames(self.Data);			% get all Data fields
   Data		= self.Data;					% copy the vector (for clarity)

   for i = 1 : length(fnames)					% for every field in Data...	
    Data.(fnames{i}) = self.Data.(fnames{i}) ( index );	% ...keep only the good data
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

 end	% methods

end	% Utils class
