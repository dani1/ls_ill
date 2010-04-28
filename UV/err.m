A = 0.01 : 0.01 : 2;
dT = 0.3;  
dA = dT ./ (100*log(10)*10.^(-A));

plot( A, 100 * dA ./ A, 'linewidth', 2, 'color', [0.5 0 1] )