classdef Experiment < dynamicprops
%============================================================================
%
% Experiment class for plotting, fitting, showing dynamic light scattering data
% taken from one Data class or from a vector of classes.
%
% Author:	Fabio Zanini
% Companies:	Institut Laue Langevin, Universität Tübingen
% Date:		12 May 2010
% Version:	2.0
% Copyright:	2009,2010 Fabio Zanini
% License:	GPLv3
%
% Syntax: dlclass = Experiment ( vector_of_points )
%============================================================================

 %============================================================================
 % PUBLIC PROPERTIES
 %============================================================================
 properties ( Access = public )

  % general parameters
  Protein
  Unit_C									% units
  Unit_Q									% units (A^{-1})
  Unit_Q2									% units (A^{-2})
  T										% temperature ( K )
  Unit_T									% units for compressibility
  Point										% data stored as fields (vectors or cell arrays)


  % N.B.: some dynamic properties could be called later on

 end	% public properties

 %============================================================================
 % DEPENDENT PROPERTIES
 %============================================================================
 properties ( GetAccess = public, Dependent, SetAccess = private )
 % These are useful for triggering events after change in a property value
 % In this class, the Data struct is changed according to changes in Q and C

  C
  Phi
  Cs
  Angle
  Q
  Q2

 end

 %============================================================================
 % PUBLIC METHODS
 %============================================================================
 methods

  %============================================================================
  % CONSTRUCTOR
  %============================================================================
  function self	= Experiment ( classin )

   if strfind(class(classin),'Sample')
    point	= [classin.Point];
   else
    point	= [classin];
   end

   % eat the common strings from the first Data class
   self.Protein		= point(1).Protein;
   self.Unit_C		= point(1).Unit_C;
   self.Unit_Q		= point(1).Unit_Q;
   self.Unit_Q2		= regexprep(self.Unit_Q,'\^\{*-1\}*','\^\{-2\}');
   self.Unit_T		= point(1).Unit_T;

   self.Point		= point;						% import the point classes as props

  end	% constructor

 end

 methods

  %============================================================================
  % GET METHODS FOR DEPENDENT PROPERTIES
  %============================================================================
  function C = get.C ( self )
   C	= unique( [ self.Point.C ] );
  end

  function Phi = get.Phi ( self )
   Phi	= unique( [ self.Point.Phi ] );
  end

  function Cs = get.Cs ( self )
   Cs	= unique( [ self.Point.Cs ] );
  end

  function Angle = get.Angle ( self )
   Angle	= unique( [ self.Point.Angle ] );
  end

  function Q = get.Q ( self )
   Q	= unique( [ self.Point.Q ] );
  end

  function T = get.T ( self )
   T	= mean( [ self.Point.T ] );
  end

 end

 % FIT METHODS
 methods

  function fit ( self, method )
   for i = 1 : length( self.Point )
    self.Point(i).fit( method );
   end
  end

  function invert_laplace ( self )
   for i = 1 : length( self.Point )
    fprintf([num2str(i) ': ']);
    self.Point(i).invert_laplace;
    fprintf('\n');
   end
  end

  function cn = coeffnames ( self, method )
   cn	= coeffnames(self.Point(1).(['Fit_' method]));
  end

  function cf = coeffvalues ( self, method )
   for i = 1 : length(self.Point)
    cf(:,i)		= coeffvalues(self.Point(i).(['Fit_' method]))';
   end
  end

 end

end
