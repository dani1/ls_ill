classdef Point < dynamicprops
% This class describes a single SLS datapoint

 properties

  Protein					% what protein?
  Salt						% what salt?
  C						% protein concentration
  Unit_C	= 'g/l'
  Cs						% salt concentration
  Unit_Cs	= 'mM'

  Angle
  Unit_Q	= 'A^{-1}'

  T
  Unit_T	= 'K'
  n

  Tau
  G
  dG

 end

 properties ( Dependent, SetAccess = private )

  Phi						% volume fraction
  Q

 end

 properties ( Hidden )

  Instrument
  C_set
  n_set

 end

 properties ( Hidden )

  Tau_raw
  G_raw
  dG_raw

 end

 % VARIOUS METHODS
 methods

  function Phi = get.Phi ( self )
   v0		= LIT.(self.Protein).v0;
   Phi		= self.C * v0;
  end

  function Q = get.Q ( self )
   Qexp 	= '4 * pi * n * sind( 0.5 * theta ) / lambda';
   Qf		= inline(Qexp,'n','theta','lambda');
   Q		= Qf( self.n, self.Angle, self.Instrument.Lambda );
  end

  function correct_G ( self )

   i	= ( self.Tau_raw < 1e-5 );			% limit time: 10 ns
   norm	= abs(mean( self.G_raw( i ) ));

   i	= ( self.Tau_raw > 1e-3	);			% shortest time: 1 mus
   i	= i & ( self.Tau_raw < 1e2	);		% longest time: 100 ms

   self.Tau	= self.Tau_raw (i);
   self.G	= self.G_raw(i)		/ norm;		% normalize
   self.dG	= self.dG_raw(i)	/ norm;

  end

 end

 % STATIC METHODS
 methods ( Static )

  fit_obj 	= fit_discrete ( t, g, dg, method, q, protein)
  [s, g, b] 	= contin  ( t, y, var, s0, s1, m, alpha, kernel)
  [ s g ]	= contin2 ( t, gt, dg, smin, smax, m, alpha, cycles )

 end

 % FIT METHODS
 methods

  function fit ( self, method )

   fit_obj	= self.fit_discrete ( self.Tau, self.G, self.dG, method, self.Q, self.Protein );

   try self.addprop(['Fit_' method]);	end
   self.(['Fit_' method])	= fit_obj;

  end

  function invert_laplace ( self )

   % PART 1: FILTER THE DATA
   ind		= ( self.Tau > 1e-3 & self.Tau < 50 & self.G > 0 );
   Tau		= self.Tau(ind);
   G		= self.G(ind);
   dG		= self.dG(ind);

   % PART 2: REDUCE THE DATA
   l		= length(Tau);
   N		= 10;
   lN		= floor(l/N);
   closeto	= 5;

   for i = 1 : N
    ind		= [ lN * (i-1) + 1 : lN * (i -1 ) + closeto ];
    t(i)	= mean( Tau(ind) );
    gt(i)	= mean( G(ind) );
    dgt(i)	= mean( dG(ind) );
   end

   y		= sqrt(gt);
   dy		= 0.5 ./ y .* dgt;

   % PART 3: PERFORM INVERSE LAPLACE TRANSFORM
   [ s, gs, bs ]= self.contin(t, y, dy, min(t), max(t), 15*N, 0.1, 0);
%   [ s, gs, bs ]= self.contin(t, y, dy, min(t), max(t), N, 0.1, 0);
   D	= 1e-6 ./ ( self.Q^2 * s );

   try   self.addprop('CONTIN'); end			% maybe it is already a property
   self.CONTIN.S	= s;
   self.CONTIN.D	= D;
   self.CONTIN.Gs	= gs;

  end

 end

end
