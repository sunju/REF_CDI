function test_bs_2D()	

	clear all; clc; 

	addpath(genpath('manopt')); 

	n = 150; 
%	x = randn(n, n);    % groundtruth 2D signal 

	x = im2double(rgb2gray(imread('mimivirus.png'))); 
	x = imresize(x, [n, n]); 
	
	
	over_fact = 6;   % over sampling factor 
	obsv_dim = round(n*over_fact); 
	bs_num = 30;     % num of beam stopping pixels in each dim 
	rmv_start = floor(obsv_dim/2) - floor(bs_num/2); 
	rmv_end = floor(obsv_dim/2) - floor(bs_num/2) + bs_num -1;
%	rmv_op = @(z) [z(1:rmv_start-1); zeros(bs_num, 1); z(rmv_end+1:end)];
	function  result = rmv_op(z)
	    result = z; 
		result(rmv_start:rmv_end, rmv_start:rmv_end) = 0; 
	end 

	y = fftshift(abs(fft2(x, obsv_dim, obsv_dim)).^2); 

	x_corr = conv2(fliplr(flipud((x))), x);   % groundtruth auto-correlation 

	aug_corr = @(z) [zeros(obsv_dim, ceil((obsv_dim - 2*n+1)/2)), [zeros(ceil((obsv_dim - 2*n+1)/2), size(z, 2)); z; zeros(floor((obsv_dim - 2*n+1)/2), size(z, 2))], zeros(obsv_dim, floor((obsv_dim - 2*n+1)/2))]; 
	ext_corr = @(v) v(ceil((obsv_dim - 2*n+1)/2)+1: ceil((obsv_dim - 2*n+1)/2) + 2*n-1, ceil((obsv_dim - 2*n+1)/2)+1: ceil((obsv_dim - 2*n+1)/2) + 2*n-1); 
	center_fft = @(z)  fftshift(fft2(ifftshift(z))); 
	center_ifft = @(z)  fftshift(ifft2(ifftshift(z))); 

	% Create the problem structure.
	manifold = euclideancomplexfactory(2*n-1, 2*n-1);
	problem.M = manifold;

	% Define the problem cost function and its Euclidean gradient.
	problem.cost  = @(v) 1.0/n^2*norm(rmv_op(y - center_fft(aug_corr(v))), 'fro').^2;
	problem.egrad = @(v) -2.0/n^2*ext_corr(center_ifft(rmv_op(y - center_fft(aug_corr(v)))));      % notice the 'e' in 'egrad' for Euclidean
	problem.ehess = @(v, w) 2.0/n^2*ext_corr(center_ifft(rmv_op(center_fft(aug_corr(w))))); 
 

	% Solve.
	options.tolgradnorm = 1e-16*n; 
    options.maxiter = 2e3; 
    options.Delta_bar = n*sqrt(problem.M.dim()); 
	[sol, sol_cost, info, options] = trustregions(problem, [], options);
	
	norm(imag(sol), 'fro')
	
	sol = real(sol); 
	
	norm(sol - 1/2*(sol + fliplr(flipud(sol))), 'fro')
	
	sol = 1/2*(sol + fliplr(flipud(sol))); 

	norm(x_corr - sol, 'fro')/norm(x_corr, 'fro') 
	 
end 
