close all; clear all; clc; 

addpath(genpath('manopt')); 

n = 10; 
over_fact = 4; 
obsv_dim = round(n*over_fact); 

x = randn(n, 1); 

y = abs(fft(x, obsv_dim)).^2; 

x_corr = conv(flip(x), x); 

pad_corr =@(v_corr) 

% Create the problem structure.
manifold = euclideancomplexfactory(2*n-1);
problem.M = manifold;
 
sub_vec = @(v) v(1:(2*n-1)); 
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(v) norm(y - fft([zeros(ceil((obsv_dim - 2*n+1)/2), 1); v; zeros(floor((obsv_dim - 2*n+1)/2), 1)]), 'fro').^2;
problem.egrad = @(v) -2*sub_vec(ifft(y - fft(v, obsv_dim)));      % notice the 'e' in 'egrad' for Euclidean
 
% Numerically check gradient consistency (optional).
checkgradient(problem);
 
% Solve.
[sol, sol_cost, info, options] = trustregions(problem, randn(2*n-1, 1));
 
 sol 
% Display some statistics.
figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');



