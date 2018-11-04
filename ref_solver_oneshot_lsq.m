function [z,err]=ref_solver_oneshot_fast(tol,maxit,solver,noise_model,img,ref,ref_type,len,centersize,alpha_min,alpha_steps,alpha_max,num_trials)

rng(1);
nimg=size(img);
nref=size(ref);
k=centersize;
alpha=1;

x = [img, ref];
n = size(x);
xpad = zeros(len);
xpad(1:n(1), 1:n(2)) = x;
f = fft2(xpad);
y = abs(f).^2; 
%% Add noise
if noise_model == 0
end
if noise_model == 1
    y=poissrnd(y);
end
if noise_model == 2
    n_photon_order = 9;
    n_photon = 1.67 * 10^n_photon_order;
    nor_fac = max(abs(f(:)));
    f = nor_fac * sqrt(  n_photon^-1 * poissrnd( n_photon/nor_fac^2 * abs(f).^2 ));
    y = f.^2;
end
%% Get autocorrelation
if solver == 0 
    r=ifft2(y);
    r=circshift(r,n(1)-1,1);
    r=circshift(r,n(2)-1,2);
    r = r(1:nimg(1),1:nimg(2));
end        
if solver == 1
    %% Remove center frequencies from data
    four_supp1 = 1+(k(1)-1)/2+1:len(1) - (k(1)-1)/2;
    four_supp2 = 1+(k(2)-1)/2+1:len(1) - (k(2)-1)/2;
    yt = y(four_supp1, four_supp2);
    %% Get autocorrelation
    vec_yt = vec(yt);
    vec_r = fft_handle_solver2d_fast(tol, maxit, vec_yt, n,len, four_supp1, four_supp2,centersize);
    r = reshape(vec_r, 2*n(1) - 1,2*n(2) - 1);
    r = r(1:nimg(1),1:nimg(2));
end
%% Run algorithm
z=img_recov(r, nimg, alpha, ref, ref_type);
%% Observed error
err=norm(img(:)-z(:))/norm(img(:));