% This script tries to apply prescriptions from
% Klotz and Walker, The Binding of Organic Ions and Proteins,
% J.Am.Chem.Soc., 1948, 70(9), pp.2935-2941
function r = ion_binding ( x, k_intr, method )

 n = 10;					% number of binding sites
% x = 1e-3;					% molar concentration of salt (M)
% k_intr = 2e6;					% decent trial constant

 if nargin < 3
  method = 'DH';
 end
 
switch method

 case 'free'
  % REACTION CONSTANTS WITHOUT INTERACTIONS (ENTROPY ONLY)
  for i = 1 : n
   k(i)	= (n - (i-1)) / i * k_intr;		% no interactions
  end
 
 case 'DH'
  % REACTION CONSTANTS WITH SCREENED COULOMB (DEBYE-HUECKEL)
  Na	= 6.02e23;				% Avogadro number
  kb	= 1.38e-23;				% Boltzmann constant ( J / K )
  T	= 300;					% temperature (Kelvin)
  z	= 1;					% charge of the ions
  e	= 1.06e-19;				% electronic charge ( C )
  eps	= 80 * 8.85e-12;			% permettivity of water ( F/m )
  b	= 30e-10;				% radius of the protein molecule ( m )
  a	= 4e-10;				% distance of minimal approach for the ion ( m )
  I	= 0.5 * ( z^2 + z ) * x;		% ionic strength of the solution
  kappa	= sqrt( (2*Na*e^2*I) / (eps*kb*T) );	% inverse Debye length
  
  %k_intr	= 2e6;
  k(1)	= k_intr;
  for i = 2 : n
   k(i) = k(i-1) * (n-i+1)/(n-i+2)*(i-1)/i*exp( - (z*e)^2/(eps*kb*T) * ( 1/b - kappa/(1+kappa*a)  ) );
  end
 
end

 % CALCULATE RATIO
 h = 0;						% initialize r
 for i = 1 : n
 
  ki	= prod(k(1:i));				% product of the first i-th reaction constants
  h	= i * ki * x^i;				% total sum
 
 end
 
 r = h / ( 1 + h );				% moles of bound ion over moles of bound protein
 
end
