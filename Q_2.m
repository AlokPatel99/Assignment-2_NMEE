%% ------------------ ESCE:543-Numerical Methods for Electrical Engineering ----------------- %%
% Assignmnet-2: Q-2: (a),(b)&(c) - Solving Potential and the Capacitance Calculations.
% The file which is used to read to data for the calculations, and it is named as input_data.txt
% The function is called from the other file provided to us named as SIMPLE2D_M.m
% The output of this file is the potential at (0.06,0,04) viz. 16th node.
% Other output is the capacitance calculation per unit lenght between the conductor and ground. 
%% ----- Code Start(Calculating the potential) ----- %%
clear
clc
%Potential = SIMPLE2D_M('xy_data.txt','sequence_data.txt','boundary_data.txt');
Potential = SIMPLE2D_M('input_data.txt');
aa = 'The Potential at the point (0.06,0.04) is: ';
bb = Potential(16,:); cc = ' V'; 

disp('-------------------------------------------');
disp('  |  Node  |     x    |    y   | Voltage|')
disp('-------------------------------------------');
disp(Potential)
disp('----------------------------------------------------------');
disp([aa num2str(bb(4)) cc])
disp('----------------------------------------------------------');

%% ----- Calculating the capacitance per unit length ----- %%

fid = fopen('input_data.txt','r');
data = textscan(fid, '%f%f%f%f',  'CollectOutput', 1);
fclose(fid);
data = data{1};
%disp(data)
Sa = [1 -0.5 -0.5;
     -0.5 0.5 0;
     -0.5 0 0.5];
Sb = [0.5  -0.5 0;
     -0.5 1 -0.5;
      0 -0.5 0.5];
E = zeros(46,1);
U = zeros(3,1);
V = 110;
vert = 35;      % 35 because in data the vertex data is stored from that position.
Energy = 0;
eps = 8.85e-12;
for i=1:46
    if mod(i,2)~=0
        S = Sa;
        %disp('Odd')
    else
        S = Sb;
        %disp('Even')
    end
    vertex = data(vert,1:3);
    %disp(vertex)
    U(1,1) = Potential(vertex(1),4);
    U(2,1) = Potential(vertex(2),4);
    U(3,1) = Potential(vertex(3),4);
    %disp(U)
    %disp(S)
    E(i) = (U'*S*U);
    %disp(E(i))
    Energy = Energy + (E(i)/2);
    vert = vert+1;
end
Ener = 4*eps*Energy;          %4 times because converted 1/4th.
C = (2*Ener)/(V^2);
ee = ' F';
dd = 'The capacitance per unit length of the system is: ';
disp([dd num2str(C) ee])
disp('-----------------------------------------------------------');

%% ----- Code End ----- %%