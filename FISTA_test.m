clear;
clc;
close all;

dim = 1000;
n   = 600;


%%%
rng(1024);
%%%

A = randn(n, dim);
x = randn(dim, 1);
x = (abs(x)>0.5) .* x;
y = A * x; 
y = y + randn(size(y)) * 0.0001;

d = 100; 


rho1 = 5;
maxIter = 300;
tol = 1e-8;

disp('Optimization')
[x_star, funcVal] = FISTA_example(A, y, rho1, maxIter, tol);
f_star = funcVal(end);

figure;
hold on;
plot(funcVal)
disp(funcVal(end))


%[x,x_star]


fprintf('x density: %.4f\n',      nnz(x)/dim)
fprintf('x_star density: %.4f\n', nnz(x_star)/dim)