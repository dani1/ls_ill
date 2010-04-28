% Example script
% It is assumed that this file is in the 'example' folder

% This is needed to execute the UV script to calculate the protein concentration
addpath('..');
addpath('../UV');

% LS paths
LSpath{1}	= 'data/LS/2010_03_27_FZ/BSA_100m_1gl_8.3mMYCl3';
LSpath{2}	= 'data/LS/2010_03_27_FZ/BSA_100m_4gl_8.3mMYCl3';
LSpath{3}	= 'data/LS/2010_03_28_FZ/BSA_100m_15gl_8.3mMYCl3';

% concentrations set in the LS instrument
c_set(1) = 1;
c_set(2) = 4;
c_set(3) = 15;

% UV paths
UVpath{1}	= 'data/UV/100328/BSA_c1_g1_8.3mMYCl3-1.jws.txt';
UVpath{2}	= 'data/UV/100328/BSA_c4_g1over4_8.3mMYCl3-1.jws.txt';
UVpath{3}	= 'data/UV/100328/BSA_c15_g1over15_8.3mMYCl3-1.jws.txt';

% pipettes used for the pre-UV dilution
UVpip{1} = { 600, {'b'}, 1 };
UVpip{2} = { [ 200 600 ], {'y' 'b'}, 1 };
UVpip{3} = { [ 70 980 ], {'y' 'b'}, 1 };

% calculate the concentrations by UV absorption at 280nm
c = c_dataset(UVpath,UVpip);

% create the Data class
for i = 1 : length(c)
  fprintf(['Example: Class n.',num2str(i),'\n']);
  data(i)	= Data(LSpath{i},'C_set',c_set(i),'C',c(i),'Repetitions',3);
end

% create the SLS and the DLS classes
sls = SLS(data);
dls = DLS(data);

clear LSpath UVpath UVpip c_set c i
