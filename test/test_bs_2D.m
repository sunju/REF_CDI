function test_bs_2D()

	clear all; clc;

	addpath(genpath('manopt'));

	% n = 150;
%	x = randn(n, n);    % groundtruth 2D signal

	% x = im2double(rgb2gray(imread('mimivirus.png')));
	% x = imresize(x, [n, n]);

	nh = 100;         % total height
	nw  = 200; 				% total width
	x = randn(nh, nw);

	mh = 1000;      % observation height
	mw = 1000; 			% observation width

	bsh  = 40; 			% beam stop height
	bsw  = 40; 			% beam stop width

	rmv_start_h = floor(mh/2) - floor(bsh/2);
	rmv_end_h = floor(mh/2) - floor(bsh/2) + bsh -1;
	rmv_start_w = floor(mw/2) - floor(bsw/2);
	rmv_end_w = floor(mw/2) - floor(bsw/2) + bsw -1;

	function  result = rmv_op(z)
	    result = z;
			result(rmv_start_h:rmv_end_h, rmv_start_w:rmv_end_w) = 0;
	end

	y = fftshift(abs(fft2(x, mh, mw)).^2);

	x_corr = conv2(fliplr(flipud((x))), x);   % groundtruth auto-correlation

	aug_corr = @(z) [zeros(mh, ceil((mw - 2*nw+1)/2)), [zeros(ceil((mh - 2*nh+1)/2), size(z, 2)); z; zeros(floor((mh - 2*nh+1)/2), size(z, 2))], zeros(mh, floor((mw - 2*nw+1)/2))];
	ext_corr = @(v) v(ceil((mh - 2*nh+1)/2)+1: ceil((mh - 2*nh+1)/2) + 2*nh-1, ceil((mw - 2*nw+1)/2)+1: ceil((mw - 2*nw+1)/2) + 2*nw-1);
	center_fft = @(z)  fftshift(fft2(ifftshift(z)));
	center_ifft = @(z)  fftshift(ifft2(ifftshift(z)));

	% Create the problem structure.
	manifold = euclideancomplexfactory(2*nh-1, 2*nw-1);
	problem.M = manifold;

 N = nw*nh;
	% Define the problem cost function and its Euclidean gradient.
	problem.cost  = @(v) 1.0/N*norm(rmv_op(y - center_fft(aug_corr(v))), 'fro').^2;
	problem.egrad = @(v) -2.0/N*ext_corr(center_ifft(rmv_op(y - center_fft(aug_corr(v)))));      % notice the 'e' in 'egrad' for Euclidean
	problem.ehess = @(v, w) 2.0/N*ext_corr(center_ifft(rmv_op(center_fft(aug_corr(w)))));


	% Solve.
	options.tolgradnorm = 1e-16*sqrt(N);
    options.maxiter = 2e3;
    options.Delta_bar = sqrt(N)*sqrt(problem.M.dim());

% initialization using zero-filled inverse Fourier transform result
		sol0 = ext_corr(center_ifft(rmv_op(y)));

	[sol, sol_cost, info, options] = trustregions(problem, sol0, options);

	norm(imag(sol), 'fro')

	sol = real(sol);

	norm(sol - 1/2*(sol + fliplr(flipud(sol))), 'fro')

	sol = 1/2*(sol + fliplr(flipud(sol)));

	norm(x_corr - sol, 'fro')/norm(x_corr, 'fro')

end
