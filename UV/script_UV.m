%==========================================================================
%
% Call the script generating the various preparations
%
%==========================================================================
dilution_table;

%==========================================================================
%
% Incluse all files into the script initial error is zero, but it is calculated below
%
%==========================================================================
% protocol with no salt
%path = '../meshdata/091127-28-Marcus/';
%protocol = { 	struct( 'file', 'BSA_1mgml_200nm-1.jws.txt', 'concentration', 1, 'error_concentration', 0, 'dilution', g(1),    'error_dilution', dg(1)) ...
%		struct( 'file', 'BSA_2mgml_200nm-1.jws.txt', 'concentration', 2, 'error_concentration', 0, 'dilution', g(2),    'error_dilution', dg(2)) ...
%		struct( 'file', 'BSA_3mgml_200nm-1.jws.txt', 'concentration', 3, 'error_concentration', 0, 'dilution', g(3),    'error_dilution', dg(3)) ...
%		struct( 'file', 'BSA_4mgml_200nm-1.jws.txt', 'concentration', 4, 'error_concentration', 0, 'dilution', g(4),    'error_dilution', dg(4)) ...
%		struct( 'file', 'BSA_5mgml_200nm-1.jws.txt', 'concentration', 5, 'error_concentration', 0, 'dilution', g(5),    'error_dilution', dg(5)) ...
%		struct( 'file', 'BSA_6mgml_200nm-1.jws.txt', 'concentration', 6, 'error_concentration', 0, 'dilution', g(6),    'error_dilution', dg(6)) ...
%		};

% protocol with 100mM Acetate buffer
path = '../../UVVis/meshdata/091210-11/';
protocol = { 	struct( 'file', 'BSA 1mgml gamma=1-1.jws.txt', 'concentration', 1, 'error_concentration', 0, 'dilution', g(11),    'error_dilution', dg(11)) ...
		struct( 'file', 'BSA 2mgml gamma=0.5-1.jws.txt', 'concentration', 2, 'error_concentration', 0, 'dilution', g(12),    'error_dilution', dg(12)) ...
		struct( 'file', 'BSA 3mgml gamma=0.3333-1.jws.txt', 'concentration', 3, 'error_concentration', 0, 'dilution', g(13),    'error_dilution', dg(13)) ...
		struct( 'file', 'BSA 4mgml gamma=0.25-1.jws.txt', 'concentration', 4, 'error_concentration', 0, 'dilution', g(3),    'error_dilution', dg(3)) ...
		struct( 'file', 'BSA 5mgml gamma=0.20-1.jws.txt', 'concentration', 5, 'error_concentration', 0, 'dilution', g(14),    'error_dilution', dg(14)) ...
		struct( 'file', 'BSA 6mgml gamma=1over6-1.jws.txt', 'concentration', 6, 'error_concentration', 0, 'dilution', g(15),    'error_dilution', dg(15)) ...
		};



dm = 0.05; % error of lab balance [mg]
dV = sqrt(2) * pipetting_error( 1000, 'blue 100-1000' ) /1000; %[ml]
n = size( protocol, 2 );
for i = 1 : n
    protocol{i}.error_concentration = sqrt( ( protocol{i}.concentration * dV ) .^2 + dm .^2 );
end

%==========================================================================
%
%  UV/VIS spectra
%
%==========================================================================

% the UV/VIS spectra measured the absorbance A = l * eps * c, where eps is
% the molar extinction coefficient which is wavelength dependent, l is path
% length and c the molar concentration, proteins strongly absorb at a
% wavelength of 279nm, the mol extinction coefficient for BSA is 
eps    = 4.42881 * 1E4; % [M^-1 * cm^-1]
length = 1;          % [cm]
mw     = 66431.5;        % molecular weight of bsa

% plot spectra of all measurements and compute molar concentration of BSA
% solution from absorbance of diluted solution at a wavelength of 279nm
concentration       = zeros( 1, n );
molar_concentration = zeros( 1, n );
error_concentration = zeros( 1, n );
error_molar_concentration = zeros( 1, n );

for i = 1 : n

    % read the data from files
    fid = fopen( [ path, protocol{i}.file ] );
    str		= '';
    wavelength	= [];
    absorbance	= [];
    while ~strcmpi(str,'XYDATA')
     str = fgetl(fid);
    end
    while ~strcmp(str,'') && ~feof(fid)
     str	= fgetl(fid);
     [ wl ab ]	= strread(str,'%f %f');
     wavelength	= [ wavelength; wl ];
     absorbance	= [ absorbance; ab ];
    end
    fclose(fid);

    concentration(i) = protocol{i}.concentration;
    error_concentration(i) = protocol{i}.error_concentration;
    
    A = absorbance(wavelength == 279 );
    
    % I checked the hardware specifications of the V630 spectrometer 
    % where they stated to have an error of 0.3 for the transmittance   T [%]
    % using the formula T = 10^-A, I derived the error for A
    dT = 0.3;
    dA = dT / (100 * log(10) * 10^(-A));
    
    % compute molar concentration using the Beer-Lambert-law
    molar_concentration(i) = A / (eps * length) * 1000 / protocol{i}.dilution; % [mM]
    error_molar_concentration(i) = A / (eps * length) * 1000 * protocol{i}.error_dilution / protocol{i}.dilution ^ 2 ...
                                   + dA / (eps * length) * 1000 / protocol{i}.dilution;
    
%    figure
%    set( gcf, 'color', 'white' )    
%    plot( wavelength, absorbance, 'linewidth', 2, 'color', [0.5, 0, 0.7] );
%    xlabel( 'wavelength [nm]', 'fontweight', 'bold', 'fontsize', 12 )
%    ylabel( 'absorbance', 'fontweight', 'bold', 'fontsize', 12 )
%    title( ['BSA ', num2str(protocol{i}.concentration),'mg/ml, (dilution factor = ',num2str(protocol{i}.dilution),')'] );
%    saveas( gcf, ['../figures/UVVis_spectrum_BSA_',num2str(protocol{i}.concentration),'_dilution_',num2str(protocol{i}.dilution),'.pdf'])
end

%==========================================================================
%
% fit data to determine specific volume of BSA
%
%==========================================================================
%
%indices = 1 : n;
%
%s = fitoptions( 'Method','NonlinearLeastSquares',...
%                'MaxFunEvals', 100000, ...
%                'Weights', 1 ./  error_molar_concentration( indices ) .^ 2, ...
%                'Lower', -inf,...
%                'Upper', inf,...
%                'Startpoint', 1 );
% 
%% we replace the keyword resolution or res with the
%% resolution at the current q value
%             
%f = fittype( '1E3 * c / (c * vbsa * 1E-3 + 1) / 66431.5', ...
%             'independent', 'c', ...
%             'coefficients', 'vbsa', ...
%             'options', s );
%                
%% fit data points
%[ cf gof ] = fit( concentration', molar_concentration', f );
%                
%% get fit result
%parameter  = coeffvalues( cf );
%confidence = confint( cf, 0.95 );
%dparameter = 0.5 * (confidence( 2, : ) - confidence( 1, : ));
%
%vbsa  =  parameter( 1 );
%dvbsa = dparameter( 1 );
%
%==========================================================================
%
% molar concentration [mM] versus concentration [mg/ml]
%
%==========================================================================
%
%figure
%hold on
%set( gcf, 'color', 'white' ) 
%set( gca, 'box', 'on', 'xscale', 'log' )
%
%xlabel( 'concentration c [mg/ml]', 'fontsize', 12 )
%ylabel( 'molar concentration c_m [mM]', 'fontsize', 12 )
%title( 'BSA molar concentration depending on concentration' )
%
%ploterr( concentration, ...
%          molar_concentration, ...
%          error_concentration, ...
%          error_molar_concentration, ...
%          'o','logx')
%
%c =  0 : range( concentration ) / 1000 : max( concentration );
%h = plot( c, 1E3 * c ./ (c * vbsa * 1E-3 + 1) / 66431.5, 'm', 'linewidth', 2 );
%
%%==========================================================================
%%
%% fit data anew including errorbars on the x-axis (balance errors),
%% propagated through the first fit, and plot everything together
%%
%%==========================================================================
%
%indices = 1 : n;
%
%new_errors = sqrt( (error_molar_concentration).^2 + ( 1E3 ./ ( (concentration .* vbsa .* 1E-3 + 1) .^2 ) ./ 66431.5 .* error_concentration ) .^2 );
%
%s = fitoptions( 'Method','NonlinearLeastSquares',...
%                'MaxFunEvals', 100000, ...
%                'Weights', 1 ./  new_errors .^ 2, ...
%                'Lower', -inf,...
%                'Upper', inf,...
%                'Startpoint', 1 );
% 
%% we replace the keyword resolution or res with the
%% resolution at the current q value
%             
%f = fittype( '1E3 * c / (c * vbsa * 1E-3 + 1) / 66431.5', ...
%             'independent', 'c', ...
%             'coefficients', 'vbsa', ...
%             'options', s );
%                
%% fit data points
%[ cf gof ] = fit( concentration', molar_concentration', f );
%                
%% get fit result
%parameter  = coeffvalues( cf );
%confidence = confint( cf, 0.95 );
%dparameter = 0.5 * (confidence( 2, : ) - confidence( 1, : ));
%
%vbsa  =  parameter( 1 );
%dvbsa = dparameter( 1 );
%Unit_vbsa = 'ml/g';
%
%
%% plot together
%errorbar ( concentration, ...
%          molar_concentration, ...
%          error_molar_concentration, ...
%          'o', 'color', 'r' )
%
%h = plot( c, 1E3 * c ./ (c * vbsa * 1E-3 + 1) / 66431.5, 'linewidth', 2, 'color', 'g' );
%
%text( 10, 3.5, '$c_m = \frac{c}{\textrm{mw} \left(1+v c\right)}$', 'interpreter','latex', 'fontsize', 18 )
%legend( 'boxoff');
%legend('hide');
%legend( h, ['specific volume v = (', num2str(vbsa), '\pm', num2str(dvbsa), ') ml/g'], 'interpreter','latex','Location','West' );
%saveas( gcf, '../figures/molecular_concentration.pdf' )

