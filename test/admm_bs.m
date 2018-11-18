
clear all; close all;

nh = 50;      % image height
nw = 50;		 % image width
mh = 6*nh;     % measurement height
mw = 6*nw;		 % measurement width

image_file = 'cat.jpg';
X0 = imresize(im2double(rgb2gray(imread(image_file))), [nh, nw]);
X0 = X0/norm(X0, 'fro');

Y  = abs(fft2(X0,mh,mw)); %abs(F'*x0*conj(G));

bsh  = 8; 			% beam stop height
bsw  = 8; 			% beam stop width

rmv_start_h = floor(mh/2) - floor(bsh/2);
rmv_end_h = floor(mh/2) - floor(bsh/2) + bsh -1;
rmv_start_w = floor(mw/2) - floor(bsw/2);
rmv_end_w = floor(mw/2) - floor(bsw/2) + bsw -1;

on_mask = ones(mh, mw);
on_mask(rmv_start_h:rmv_end_h, rmv_start_w:rmv_end_w) = 0;

% try to initialize with the right magnitude with random phases
%Z = randn(m1, m2); %initialization for Z and W
W = ifft2(Y.*on_mask);
Z = rand(mh, mw);
lambda = randn(mh, mw);

eq_gap = @(Z, W) norm(Z - W, 'fro');

Ps = @(z) [z(1:nh, 1:nw), zeros(nh, mw-nw);  zeros(mh-nh, mw)];

Pm = @(z)    ifft2(Y.*ifftshift(fftshift(exp(1i*angle(fft2(z)))).*on_mask) + ifftshift(fftshift(fft2(z)).*(1-on_mask)));

% Pm = @(z)  ifft2(Y.*exp(1i*angle(fft2(z)) ) );

GAP_TOL = 1e-7;
rho = 20;         % this might need to be tuned larger to see convergence
iter = 1;

while eq_gap(Z, W) >= GAP_TOL

	Z = Ps(W - lambda/rho);
	W = Pm(Z + lambda/rho);

	% update the dual variable
	lambda=lambda + rho*(Z-W);

	disp(['Iter = ', num2str(iter), ', Gap = ' num2str(eq_gap(Z, W))]);
	iter = iter + 1;
end


Z=real(Z);
W=real(W);

%save('fpr_2d_admm_rand');
figure;
subplot(1, 2, 1);
imagesc(X0);
subplot(1, 2, 2);
imagesc(Z(1:nh, 1:nw));

norm(abs(fft2(X0, mh, mw)) - abs(fft2(Z/norm(Z, 'fro'), mh, mw)), 'fro')/norm(abs(fft2(X0, mh, mw)), 'fro')
