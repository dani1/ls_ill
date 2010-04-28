%==========================================================================
%
% here come the various sample preparations
% each solution was diluted with water, the dilution factor is computed as
% follows: x ml BSA solution of concentration c, was mixed with y ml of H2O
% hence the dilution factor g = x / (x + y)
% all quantities are in mul
%
%==========================================================================
% 1.
s1{1}	= 200;
sn1{1}	= { 'y' };
num1{1}	= 1;

% 2.
s1{2}	= [ 150 150 150 150 ];
sn1{2}	= { 'y' 'y' 'y' 'y' };
num1{2}	= 2;

% 3.
s1{3}	= [ 200 600 ];
sn1{3}	= { 'y' 'b' };
num1{3}	= 1;

% 4.
s1{4}	= [ 150 150 900 ];
sn1{4}	= { 'y' 'y' 'b' };
num1{4}	= 2;

% 5.
s1{5}	= [ 200 800 ];
sn1{5}	= { 'y' 'b' };
num1{5}	= 1;

% 6.
s1{6}	= [ 200 1000 ];
sn1{6}	= { 'y' 'b'  };
num1{6}	= 1;

% 7.
s1{7}	= [ 212 788 ];
sn1{7}	= { 'b' 'b' };
num1{7}	= 1;
s2{7}	= [ 100 900 ];
sn2{7}	= { 'y' 'b' };
num2{7}	= 1;

% 8.
s1{8}	= [ 100 900 ];
sn1{8}	= { 'y' 'b' };
num1{8}	= 1;
s2{8}	= [ 100 900 ];
sn2{8}	= { 'y' 'b' };
num2{8}	= 1;

% 9.
s1{9}	= [ 100 900 ];
sn1{9}	= { 'y' 'b' };
num1{9}	= 1;
s2{9}	= [ 100 900 ];
sn2{9}	= { 'y' 'b' };
num2{9}	= 1;
s3{9}	= [ 180 420 ];
sn3{9}	= { 'y' 'b' };
num3{9}	= 1;

% 10.
s1{10}	= [ 100 900 ];
sn1{10}	= { 'y' 'b' };
num1{10}	= 1;
s2{10}	= [ 100 900 ];
sn2{10}	= { 'y' 'b' };
num2{10}	= 1;
s3{10}	= [ 500 500 ];
sn3{10}	= { 'b' 'b' };
num3{10}= 1;

% 11.
s1{11}	= 650;
sn1{11}	= { 'b' };
num1{11}= 1;

% 12.
s1{12}	= [ 500 500 ];
sn1{12}	= { 'b' 'b' };
num1{12}= 1;

% 13.
s1{13}	= [ 150 150 600 ];
sn1{13}	= { 'y' 'y' 'b' };
num1{13}= 2;

% 14.
s1{14}	= [ 150 600 ];
sn1{14}	= { 'y' 'b' };
num1{14}= 1;

% 15.
s1{15}	= [ 150 750 ];
sn1{15}	= { 'y' 'b' };
num1{15}= 1;

% calculate the dilutions
for i = 1 : length(s1)
 if	isempty(num2{i})
  [ g(i) dg(i) ] = dilution(s1{i},sn1{i},num1{i});
 elseif	isempty(num3{i})
  [ g(i) dg(i) ] = dilution(s1{i},sn1{i},num1{i},s2{i},sn2{i},num2{i});
 else
  [ g(i) dg(i) ] = dilution(s1{i},sn1{i},num1{i},s2{i},sn2{i},num2{i},s3{i},sn3{i},num3{i});
 end
end

