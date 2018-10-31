function [alpha_vec,err_vec]=ref_solver(solver,noise_model,img,ref,ref_type,len,centersize,alpha_min,alpha_steps,alpha_max,num_trials)

rng(1);
nimg=size(img);
nref=size(ref);
k=centersize;

i=0;
ref_unscaled=ref;
for alpha=logspace(alpha_min,alpha_max,alpha_steps)%alpha_min:alpha_steps:alpha_max
    i=i+1;
    alpha_vec(i)=alpha;
    ref=alpha*ref_unscaled;
    %%%
    x = [img, ref];
    n = size(x);
    xpad = zeros(len);
    xpad(1:n(1), 1:n(2)) = x;
    f = fft2(xpad);
    y = abs(f).^2; 
    err_vec(i)=0;
    for j=1:num_trials
        [i,j]
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
            four_supp1 = ceil(k(1)/2) + 1 : len(1) - floor(k(1)/2);
            four_supp2 = ceil(k(2)/2) + 1 : len(1) - floor(k(2)/2);
            yt = y(four_supp1, four_supp2);
            %% Get autocorrelation
            vec_yt = vec(yt);
            vec_r = fft_handle_solver2d(vec_yt, n,len, four_supp1, four_supp2);
            r = reshape(vec_r, 2*n(1) - 1,2*n(2) - 1);
            r = r(1:nimg(1),1:nimg(2));
        end
        %% Run algorithm
        z=img_recov(r, nimg, alpha, ref, ref_type);
        %% Observed error
        err=norm(img(:)-z(:))/norm(img(:));
        err_vec(i)=err_vec(i)+err;
    end
    err_vec(i)=err_vec(i)/num_trials;
end