function example()
    %example file to understand how to load data with the program.

 main_dir = '~/Documents/tesi/Matlab/ls-ill/'; % !!! INSERT HERE THE PATH OF ls-ill  !!!
 addpath('../') %access parent folder with all class-definitions
 addpath(main_dir)
 
 
 Name		= 'dY8';

 Protein	= 'BSA';
 Salt		= 'YCl3';
 cs		= 5;		% mM

 n		= Index_Refraction.H2O;
 n_set		= n;

T		= -Constants.T0 + 22;
 % LS paths 

 DLSpath{1}	= [main_dir '/' 'example/example-data/LS/2010_03_27_FZ/BSA_100m_1gl_8.3mMYCl3'];
 DLSpath{3}	= [main_dir '/' 'example/example-data/LS/2010_03_27_FZ/BSA_100m_4gl_8.3mMYCl3'];
 DLSpath{2}	= [main_dir '/' 'example/example-data/LS/2010_03_28_FZ/BSA_100m_15gl_8.3mMYCl3'];
 disp(DLSpath{1})
 % concentrations set in the LS instrument
 c_set(1) = 1;
 c_set(2) = 4;
 c_set(3) = 15;
 
 % UV paths
 UVpath{1}	= [main_dir '/' 'example/example-data/UV/100328/BSA_c1_g1_8.3mMYCl3-1.jws.txt'] ;
 UVpath{2}	= [main_dir '/' 'example/example-data/UV/100328/BSA_c4_g1over4_8.3mMYCl3-1.jws.txt'] ;
 UVpath{3}	= [main_dir '/' 'example/example-data/UV/100328/BSA_c15_g1over15_8.3mMYCl3-1.jws.txt'] ;
 
 % pipettes used for the pre-UV dilution
UVpip{1} = { 600, {'b'}, 1 };
UVpip{2} = { [ 200 600 ], {'y' 'b'}, 1 };
UVpip{3} = { [ 70 980 ], {'y' 'b'}, 1 };


 % calculate the concentrations
 
if exist([main_dir '/../UV']);%check if directory with functions extracting c from
                  %UV data exists. Only useful for example file purpose
   addpath([main_dir '/../UV'])
   c = c_dataset(UVpath,UVpip) %use c_dataset in '../UV' to calculate c from UVpath
else
   disp('c = c_set') %if directory not there set c = c_set
   c = c_set 
end


 % call the class
 for i = 1 : length(c_set)
  sample_dls(i)	= DLS.Sample(	'Path',		DLSpath{i},	...
				'Protein',	Protein,	...
				'Salt',		Salt,		...
				'T',		T,		...
				'n',		n,		...
				'n_set',	n_set,		...
				'C',		c(i),		...
				'C_set',	c_set(i),	...
				'Cs',		cs,		...
				'Instrument',	'ALV'		);
 end
 
 assignin('caller',['d' Name],sample_dls);
 
 
% SLS files saved as _bak.txt files copying the result table of ALV instrument.
 SLSpath{1} = [DLSpath{1} '_bak.txt'];
 SLSpath{2} = [DLSpath{2} '_bak.txt'];
 SLSpath{3} = [DLSpath{3} '_bak.txt'];
 
 for i = 1 : length(c_set)
     
     sample_sls(i) = SLS.Sample('Path',		SLSpath{i},	...
				'Protein',	Protein,	...
				'Salt',		Salt,		...
				'T',		T,		...
				'n',		n,		...
				'n_set',	n_set,		...
				'C',		c(i),		...
				'C_set',	c_set(i),	...
				'Cs',		cs,		...
				'Instrument',	'ALV', ...	
                'dndc_set', 0.175 , ...
                'dndc',     0.175);
 end
  fit_example(sample_dls)
 assignin('caller',['s' Name],sample_sls);
end % end of example function, how to load data from DLS and SLS files


function fit_example(sample_dls)
	for i = 1 : length(sample_dls)
		sample_dls(i).fit('Single') % other options instead of Single are (see Point): Double, Cumulants
	end
	disp('example fit result Single of on data point, other results accessible via dY8(i).Point(j).Fit_Single')
	sample_dls(1).Point(1).Fit_Single
end
