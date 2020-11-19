%% ------------------ ESCE:543-Numerical Methods for Electrical Engineering ----------------- %%
% Assignmnet-2: Q-3: Conjugate Method to solve for the potential at points.

%% ----- Code Start ----- %%
clear
clc
%{
b = [0 0 0 0 0 0 0 0 0 0 0 0 -110 -110 -110 0 -110 0 -110]';
A = [-4 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
     1 -4 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
     0 1 -4 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
     0 0 1 -4 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
     0 0 0 2 -4 1 0 0 0 1 0 0 0 0 0 0 0 0 0;
     1 0 0 0 0 -4 1 0 0 0 1 0 0 0 0 0 0 0 0;
     0 1 0 0 0 1 -4 1 0 0 0 1 0 0 0 0 0 0 0;
     0 0 1 0 0 0 1 -4 1 0 0 0 1 0 0 0 0 0 0;
     0 0 0 1 0 0 0 1 -4 1 0 0 0 1 0 0 0 0 0;
     0 0 0 0 1 0 0 0 2 -4 0 0 0 0 1 0 0 0 0;
     0 0 0 0 0 1 0 0 0 0 -4 1 0 0 0 1 0 0 0;
     0 0 0 0 0 0 1 0 0 0 1 -4 1 0 0 0 1 0 0;
     0 0 0 0 0 0 0 1 0 0 0 1 -4 1 0 0 0 0 0;
     0 0 0 0 0 0 0 0 1 0 0 0 1 -4 1 0 0 0 0;
     0 0 0 0 0 0 0 0 0 1 0 0 0 2 -4 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 -4 1 1 0;
     0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 -4 0 1;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 -4 1;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 1 -4];
 n = 19;
%}
%{
% Below is the data for the new b and A with A as PD.
 b = [0;0;0;0;0;0;0;0;0;0;0;0;110;110;110;0;110;0;110];
 A = [4,-1,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0;
        -1,4,-1,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0;
        0,-1,4,-1,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0;
        0,0,-1,4,-1,0,0,0,-1,0,0,0,0,0,0,0,0,0,0;
        0,0,0,-2,4,0,0,0,0,-1,0,0,0,0,0,0,0,0,0;
        -1,0,0,0,0,4,-1,0,0,0,-1,0,0,0,0,0,0,0,0;
        0,-1,0,0,0,-1,4,-1,0,0,0,-1,0,0,0,0,0,0,0;
        0,0,-1,0,0,0,-1,4,-1,0,0,0,-1,0,0,0,0,0,0;
        0,0,0,-1,0,0,0,-1,4,-1,0,0,0,-1,0,0,0,0,0;
        0,0,0,0,-1,0,0,0,-2,4,0,0,0,0,-1,0,0,0,0;
        0,0,0,0,0,-1,0,0,0,0,4,-1,0,0,0,-1,0,0,0;
        0,0,0,0,0,0,-1,0,0,0,-1,4,-1,0,0,0,-1,0,0;
        0,0,0,0,0,0,0,-1,0,0,0,-1,4,-1,0,0,0,0,0;
        0,0,0,0,0,0,0,0,-1,0,0,0,-1,4,-1,0,0,0,0;
        0,0,0,0,0,0,0,0,0,-1,0,0,0,-2,4,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,4,-1,-1,0;
        0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,-1,4,0,-1;
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,0,4,-1;
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-1,4];
n = 19;
%}

% New S and c are as below, by assigning value at A and b respectively:
 b = [0;0;0;0;0;0;0;-110;-110;-110;0;-220;330;110;330;-110;220;-110;330];
 A = [18,-8,1,0,0,-8,2,0,0,0,1,0,0,0,0,0,0,0,0;
     -8,19,-8,1,0,2,-8,2,0,0,0,1,0,0,0,0,0,0,0;
     1,-8,19,-8,1,0,2,-8,2,0,0,0,1,0,0,0,0,0,0;
     0,1,-8,22,-12,0,0,2,-8,3,0,0,0,1,0,0,0,0,0;
     0,0,1,-12,18,0,0,0,3,-8,0,0,0,0,1,0,0,0,0;
     -8,2,0,0,0,19,-8,1,0,0,-8,2,0,0,0,1,0,0,0;
     2,-8,2,0,0,-8,20,-8,1,0,2,-8,2,0,0,0,1,0,0;
     0,2,-8,2,0,1,-8,20,-8,1,0,2,-8,2,0,0,0,0,0;
     0,0,2,-8,3,0,1,-8,23,-12,0,0,2,-8,3,0,0,0,0;
     0,0,0,3,-8,0,0,1,-12,19,0,0,0,3,-8,0,0,0,0;
     1,0,0,0,0,-8,2,0,0,0,19,-8,1,0,0,-8,2,1,0;
     0,1,0,0,0,2,-8,2,0,0,-8,20,-8,1,0,2,-8,0,1;
     0,0,1,0,0,0,2,-8,2,0,1,-8,19,-8,1,0,1,0,0;
     0,0,0,1,0,0,0,2,-8,3,0,1,-8,22,-12,0,0,0,0;
     0,0,0,0,1,0,0,0,3,-8,0,0,1,-12,18,0,0,0,0;
     0,0,0,0,0,1,0,0,0,0,-8,2,0,0,0,22,-8,-12,3;
     0,0,0,0,0,0,1,0,0,0,2,-8,1,0,0,-8,22,3,-12;
     0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-12,3,18,-8;
     0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,3,-12,-8,18];
n = 19;
%}
%% --- Long Method --- %%
%{
for i=1:n
    xk(i,1) = 0;
end
xk1 = xk;
rk = xk;
rk1 = xk;
var2 = xk;
var1 = xk;
var3 = xk;
var4 = xk;
var5 = xk;
pk1 = xk;
for i=1:n
    for j=1:n
        var1(i,1) = var1(i,1) + (A(i,j)*xk(j,1));
    end
end
for i=1:n
    rk(i,1) = b(i,1) - var1(i,1);
end

pk = rk;
eps = 1e-5;
N = 0;
norm_inf = [];
iteration = [];
norm_2 = [];

%while true
for k=1:19
    num1=0;
    deno1=0;
    for i=1:n
        num1 = num1 + (pk(i,1)*rk(i,1));
    end
    for i=1:n
        for j=1:n
            var2(i,1) = var2(i,1) + (A(i,j)*pk(j,1));
        end
    end
    for i=1:n
        deno1 = deno1 + (pk(i,1)*var2(i,1));
    end
    alpha = num1/deno1;
    for i=1:n
        xk1(i,1) = xk(i,1) + (alpha*pk(i,1));
    end
    for i=1:n
        for j=1:n
            var3(i,1) = var3(i,1) + (A(i,j)*xk1(j,1));
        end
    end
    for i=1:n
        rk1(i,1) = b(i,1) - var3(i,1);
    end
    num2 = 0;
    deno2 = 0;
    for i=1:n
        for j=1:n
            var4(i,1) = var4(i,1) + (A(i,j)*rk1(j,1));
        end
    end
    for i=1:n
        num2 = num2 + (pk(i,1)*var4(i,1));
    end
    for i=1:n
        var5(i,1) = var5(i,1) + (A(i,j)*pk(j,1));
    end
    for i=1:n
        deno2 = deno2 + (pk(i,1)*var5(i,1));
    end
    beta = (-1)*(num2/deno2);
    for i=1:n
        pk1(i,1) = rk1(i,1) + (beta*pk(i,1));
    end

    limit = max(abs(rk1));
    two_norm = 0;
    for i=1:n
        two_norm = two_norm + (rk1(i,1)*rk1(i,1));
    end
    if two_norm<eps
        break
    end
    norm_inf = [norm_inf; limit];
    norm_2 = [norm_2; sqrt(two_norm)];

    rk = rk1;           
    pk = pk1;
    xk = xk1;
    N = N + 1;
    iteration = [iteration; N];
end
%}
%% ----- Direct Method ----- %%

xk = zeros(n,1);
rk = b - A*xk;
pk = rk;
eps = 1e-5;
N = 0;
norm_inf = [];
iteration = [];
norm_2 = [];

while true

    num1 = pk'*rk;
    deno1 = pk'*A*pk;
    alpha = num1/deno1;
    xk1 = xk + alpha*pk;
    rk1 = b - A*xk1;
    num2 = pk'*A*rk1;
    deno2 = pk'*A*pk;
    beta = (-1)*(num2/deno2);
    pk1 = rk1 + beta*pk;
    
    limit = max(abs(rk1));
    two_norm = 0;
    for i=1:n
        two_norm = two_norm + (rk1(i,1)*rk1(i,1));
    end
    if two_norm<eps
        break
    end
    norm_inf = [norm_inf; limit];
    norm_2 = [norm_2; sqrt(two_norm)];

    rk = rk1;           
    pk = pk1;
    xk = xk1;
    N = N + 1;
    iteration = [iteration; N];
end

%% --- Plot --- %%
%{
f1 = figure;
f2 = figure;

figure(f1);
plot(iteration,norm_inf)
title('Plot of infinite norm residual vs Number of iterations')
xlabel('Number of iterations')
ylabel('Infinite norm')
hold on

figure(f2);
plot(iteration,norm_2)
title('Plot of two norm residual vs Number of iterations')
xlabel('Number of iterations')
ylabel('Two Norm')
%}
aa = 'The Potential at the point (0.06,0.04) is: ';
bb = xk1; cc = ' V';
nn = [1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19];
s = ['  ';'  ';'  ';'  ';'  ';'  ';'  ';'  ';'  ';
     '  ';'  ';'  ';'  ';'  ';'  ';'  ';'  ';'  ';'  '];
disp('-------------------------------------------');
disp('  |  Node  | Voltage|')
disp('-------------------------------------------');
disp([s s num2str(nn) s s s num2str(xk1)])
disp('-----------------------------------------------------');
disp([aa num2str(xk1(8)) cc])
disp('-----------------------------------------------------');

%% --- Code Ends --- %%