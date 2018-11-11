clear all; clc; 

addpath(genpath('manopt')); 

n = 100; 
over_fact = 4; 
obsv_dim = round(n*over_fact); 

x = randn(n, 1); 

y = fftshift(abs(fft(x, obsv_dim)).^2); 

x_corr = conv(flip(x), x); 

aug_corr = @(z) [zeros(ceil((obsv_dim - 2*n+1)/2), 1); z; zeros(floor((obsv_dim - 2*n+1)/2), 1)]; 
ext_corr = @(v) v(ceil((obsv_dim - 2*n+1)/2)+1 : ceil((obsv_dim - 2*n+1)/2) + 2*n-1); 
center_fft = @(z)  fftshift(fft(ifftshift(z))); 
center_ifft = @(z)  fftshift(ifft(ifftshift(z))); 

% Create the problem structure.
manifold = euclideancomplexfactory(2*n-1);
problem.M = manifold;

% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(v) norm(y - center_fft(aug_corr(v)), 'fro').^2;
problem.egrad = @(v) -2*ext_corr(center_ifft(y - center_fft(aug_corr(v))));      % notice the 'e' in 'egrad' for Euclidean

% Solve.
%[sol, sol_cost, info, options] = trustregions(problem);
% 
[sol, sol_cost, info, options] = rlbfgs(problem); 
 

