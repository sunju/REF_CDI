
clear all; close all;

n1 = 50;     % image height
n2 = 50;		 % image width
m1 = 2 * n1;     % measurement height
m2 = 2 * n2;		 % measurement width

image_file = 'cat.jpg';
X0 = imresize(im2double(rgb2gray(imread(image_file))), [n1, n2]);
X0 = X0/norm(X0, 'fro');

Y  = abs(fft2(X0,m1,m2)); %abs(F'*x0*conj(G));


% try to initialize with the right magnitude with random phases
%Z = randn(m1, m2); %initialization for Z and W
Z = ifft2(Y);
W = rand(m1, m2);
lambda = randn(m1, m2);

eq_gap = @(Z, W) norm(Z - W, 'fro');

Ps = @(z) [z(1:n1, 1:n2),zeros(n1, m2-n2);  zeros(m1-n1, m2)];
Pm = @(z) ifft2(Y.*exp(1i*angle(fft2(z)) ) );

GAP_TOL = 1e-7;
rho = 5;         % this might need to be tuned larger to see convergence 
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
imagesc(Z(1:n1, 1:n2));

norm(abs(fft2(X0, m1, m2)) - abs(fft2(Z/norm(Z, 'fro'), m1, m2)), 'fro')/norm(abs(fft2(X0, m1, m2)), 'fro')
