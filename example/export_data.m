% for this script to work, the actual path has to be in the folder containing this file/
% this folder has to be added to the matlab search path
% check to modify the variably main_dir in example.m before running this file
example
dls_sample = dBSAwY8p3(1);
dls_sample.fit('DoubleBKG');
[A1 dA1]=get_fit(dls_sample,'DoubleBKG','A1');
[A2 dA2]=get_fit(dls_sample,'DoubleBKG','A2');
[G2 dG2]=get_fit(dls_sample,'DoubleBKG','Gamma2');
[G1 dG1]=get_fit(dls_sample,'DoubleBKG','Gamma1');
q=[dls_sample.Qv];
q2=q.^2;
dat=[q q2 A1 dA1 G1 dG1 A2 dA2 G2 dG2]
fnam='outData.txt'; % use your file!
str='q q2 A1 dA1 G1 dG1 A2 dA2 G2 dG2';
dlmwrite(fnam,str,'delimiter','');
dlmwrite(fnam,dat,'delimiter',' ', 'precision','%.3e','-append');
