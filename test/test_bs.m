clear all; clc; 

addpath(genpath('manopt')); 

n = 500; 
over_fact = 4; 
obsv_dim = round(n*over_fact); 
bs_num = 100;     % always assumed to be even 
rmv_start = floor(obsv_dim/2) - floor(bs_num/2); 
rmv_end = floor(obsv_dim/2) - floor(bs_num/2) + bs_num -1; 
rmv_op = @(z) [z(1:rmv_start-1); zeros(bs_num, 1); z(rmv_end+1:end)]; 

x = randn(n, 1); 

y = fftshift(abs(fft(x, obsv_dim)).^2); 

x_corr = conv(flip(x), x); 

aug_corr = @(z) [zeros(ceil((obsv_dim - 2*n+1)/2), 1); z; zeros(floor((obsv_dim - 2*n+1)/2), 1)]; 
ext_corr = @(v) v(ceil((obsv_dim - 2*n+1)/2)+1 : ceil((obsv_dim - 2*n+1)/2) + 2*n-1); 
center_fft = @(z)  fftshift(fft(ifftshift(z))); 
center_ifft = @(z)  fftshift(ifft(ifftshift(z))); 

% Create the problem structure.
manifold = euclideancomplexfactory(2*n-1);
%manifold = euclideanfactory(2*n-1);
problem.M = manifold;

% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(v) norm(rmv_op(y - center_fft(aug_corr(v))), 'fro').^2;
problem.egrad = @(v) -2*ext_corr(center_ifft(rmv_op(y - center_fft(aug_corr(v)))));      % notice the 'e' in 'egrad' for Euclidean
problem.ehess = @(v, w) 2*ext_corr(center_ifft(rmv_op(center_fft(aug_corr(w))))); 
 
 
% Solve.

options.tolgradnorm = 1e-14*sqrt(n); 
options.maxiter = 2e3; 
[sol, sol_cost, info, options] = trustregions(problem, [], options);
 
%[sol, sol_cost, info, options] = rlbfgs(problem, [], options); 

norm(x_corr - sol)/norm(x_corr)

