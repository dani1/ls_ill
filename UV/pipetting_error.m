% given volume [mul] the function interpolates the pipetting error dvolume [mul] 
% for eppendorf pipettes, errors for different pipette tips can be found in 
% 'Research Pipette Instruction Manual.pdf'
function dvolume = pipetting_error( volume, tip )

    if ( strcmp( tip, 'yellow') || strcmp( tip, 'y' ) || strcmp ( tip, 'yellow 20-200') )
        V = [ 20, 0.5; ...
             100, 1.0; ...
             200, 1.2 ];
    elseif ( strcmp( tip, 'blue') || strcmp( tip, 'b' ) || strcmp( tip, 'blue 100-1000') )
        V = [ 100, 3; ...
              500, 1; ...
              1000, 0.6 ];
    else
     error('Pipette type not recognized!');
    end
     
    dvolume = interp1( V(:,1), V(:,2), volume );
end
 
