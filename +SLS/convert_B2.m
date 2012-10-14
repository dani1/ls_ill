function B2r = convert_B2 ( B2, M )

    if nargin == 1
        M = LIT.BSA.M;
        disp(['I assume the literature mass of BSA = ' num2str(M) ' Da']);
    end
    Na	= Constants.Na;
    B2r	= 1e-3 * B2 .* M.^2 / Na;		% m^3, starting from mol*l / g^2
end
