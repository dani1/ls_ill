function B2 = B2_Donnan ( cp, cs, z )
% This function calculates the Donnan contribution to the
% second virial coefficient of protein-salt-water solutions
% according to:
%
% Asthagiri et al., A Consistent Experimental and Modeling
% Approach to Light-Scattering Studies of Protein-Protein
% Interactions in Solution, Biophys. J., 88, 3300 (2005)
%
% The formula is the following:
%
B2f	= inline('z.^2 ./ ( 2 .* cs + z .* cp)','cp','cs','z');
% where cp and cs are the molar protein and salt concentrations,
% expressed in mol/l, and z is the association number which
% assures electroneutrality of the chosen components.
% It should be valid to assume z = |protein charge|

B2	= B2f(cp,cs,z);

% To get B2 in usual units, i.e. l*mol/g^2, one has to divide
% for the squared molar mass of the protein

end
