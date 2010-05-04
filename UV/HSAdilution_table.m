%==========================================================================
%
% here come the various sample preparations
% each solution was diluted with water, the dilution factor is computed as
% follows: x ml BSA solution of concentration c, was mixed with y ml of H2O
% hence the dilution factor g = x / (x + y)
% all quantities are in mul
%
%==========================================================================


%protocol for the dilution procedure of the samples: Which pipette zellow
%or blue?, Which volume Protein?, which volume diluted sample?, How many times was the sample diluted?
%protocol = {
%    struct('Pipette',[y ,b ], 'volumeProtein',[20, 30,50,100,200,300,350,370,390,400,500,], 'volumeDilution', [400,410,420,450,460], 'dilutionSteps',[1,2,3,4])
%};
protocol = [...
    struct('PipetteForProtein', ['y'], 'volumeProtein',[70] ,'PipetteForWater',['b'], 'volumeWater', [350]),...
    struct('PipetteForProtein', ['y'], 'volumeProtein',[70] ,'PipetteForWater',['b'], 'volumeWater', [350]),...
    struct('PipetteForProtein', ['y'], 'volumeProtein',[50] ,'PipetteForWater',['b'], 'volumeWater', [350]),...
    struct('PipetteForProtein', ['y'], 'volumeProtein',[30] ,'PipetteForWater',['b'], 'volumeWater', [390]),...
    struct('PipetteForProtein', ['y'], 'volumeProtein',[30] ,'PipetteForWater',['b'], 'volumeWater', [390]),...
    struct('PipetteForProtein', ['y'], 'volumeProtein',[20] ,'PipetteForWater',['b'], 'volumeWater', [380]),...    
    ];

%actually the pipette for 20\mu l was the white 2-20mu l pipette

for i = 1 : length(protocol)
   [ g(i) dg(i) ] = HSAdilution(protocol(i))
end





