function eta = ViscosityH2O( T, cs )
    % T in units of K
    % cs in units of [M] = [mol/l]
    %==========================================================================
    %
    % dependence of viscosity on temperature was determined by Cho et al. (1999)
    % "Thermal offset viscosities of liquid H2O, D2O and T2O", J.Phys. Chem. B
    % 103(11):1991-1994
    %
    % eta0(T) is valid from 280k up to 400K
    %
    %==========================================================================

    C  = 802.25336;
    a  = 3.4741 * 10^(-3);
    b  = -1.7413 * 10^(-5);
    c  = 2.7719 * 10^(-8);
    g  = 1.53026;
    T0 = 225.334;
    dT = T - T0;

    eta_0 = C * ( dT + a * dT.^2 + b * dT.^3 + c * dT.^4 ) .^ (-g);

    %==========================================================================
    %
    % correction for low salt concentrations (<0.1M) according to Desnoyers and
    % Perron (1972): "The viscosity of aqueous solutions of alkali and tetra-
    % -alkylammonium halides at 25C", Jornal of Solution Chemistry, 1(3).199-212 
    %
    %==========================================================================

    A = 0.006;
    B = 0.080;
    D = 0.007;

    eta_r = 1 + A * cs.^(1/2) + B * cs + D * cs.^2;

    eta = eta_0 .* eta_r * 1e-3; % [Pa s]
end
