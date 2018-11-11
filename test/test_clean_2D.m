clear all; clc; 

addpath(genpath('manopt')); 

n = 250; 
over_fact = 4; 
obsv_dim = round(n*over_fact); 

x = randn(n, n); 

y = fftshift(abs(fft2(x, obsv_dim, obsv_dim)).^2); 

x_corr = conv2(fliplr(flipud((x))), x); 

aug_corr = @(z) [zeros(obsv_dim, ceil((obsv_dim - 2*n+1)/2)), [zeros(ceil((obsv_dim - 2*n+1)/2), size(z, 2)); z; zeros(floor((obsv_dim - 2*n+1)/2), size(z, 2))], zeros(obsv_dim, floor((obsv_dim - 2*n+1)/2))]; 
ext_corr = @(v) v(ceil((obsv_dim - 2*n+1)/2)+1: ceil((obsv_dim - 2*n+1)/2) + 2*n-1, ceil((obsv_dim - 2*n+1)/2)+1: ceil((obsv_dim - 2*n+1)/2) + 2*n-1); 
center_fft = @(z)  fftshift(fft2(ifftshift(z))); 
center_ifft = @(z)  fftshift(ifft2(ifftshift(z))); 

% Create the problem structure.
manifold = euclideancomplexfactory(2*n-1, 2*n-1);
problem.M = manifold;

% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(v) norm(y - center_fft(aug_corr(v)), 'fro').^2;
problem.egrad = @(v) -2*ext_corr(center_ifft(y - center_fft(aug_corr(v))));      % notice the 'e' in 'egrad' for Euclidean

% Solve.
x0 = zeros(2*n-1, 2*n-1); 
[sol, sol_cost, info, options] = trustregions(problem, x0);
% 
%[sol, sol_cost, info, options] = rlbfgs(problem); 

norm(x_corr - sol, 'fro')/norm(x_corr, 'fro')
 
