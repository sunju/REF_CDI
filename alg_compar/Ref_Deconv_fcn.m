% Function implementing the Referenced Deconvolution algorithm. Input
% variables are:
% n - dimensions of n x n specimen and n x n reference
% m - dimensions of m x m squared Fourier transform magnitude data
% img - ground-truth specimen image
% ref_type - reference type: 'b' - block reference, 's' - slit reference,
%   'p' - pinhole reference, 'a' - arbitrary reference
% ref - reference R
% num_trials - number of trials to repeat (each trial generates a
%   different random instance of Poisson shot noise). Final results are
%   averaged over all trials
% photon_param_range - range of photon per pixel values to test
% photon_param_rec - single photon per pixel value for which the recovered
% image is saved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err_vec,exp_err_vec,img_rec]=Ref_Deconv(n,m,img,ref_type,ref,num_trials,photon_param_range,photon_param_rec)

rng(1)

%% DFT matrices and partial DFT matrices
F1=dftmtx(m);
F1c=F1(:,[end-(n-1)+1:end,1:n]);
F1cc=F1(:,1:n);
F2=dftmtx(m);
F2c=F2(:,[end-(2*n-1)+1:end,1:2*n]);
F2cc=F2(:,1:n);

%% Reference choice

if ref_type=='b';
    ref=ones(n);
    % Scaling factor
    S_dec = 1/(m^2) + 2*((n-1)/(m^2))*(1-cos(2*pi*[0:m-1]/m));
    S=S_dec'*S_dec;
end

if ref_type=='s';
    ref=zeros(n); ref(:,end)=ones(n,1);
    % Scaling factor
    S_dec = 1/(m^2) + 2*((n-1)/(m^2))*(1-cos(2*pi*[0:m-1]/m));
    S_flat = (n/(m^2))*ones(1,m);
    S=S_dec'*S_flat;
end

if ref_type=='p';
    ref=zeros(n); ref(n,n)=1;
    % Scaling factor
    S_flat = (n/(m^2))*ones(1,m);
    S=S_flat'*S_flat;
end

if ref_type=='a'
    MR=ref2mtrx(ref);
    TR=(1/(m^2))*inv(MR)*(kron(F2cc',F1cc'));
    S=reshape(diag(TR'*TR),m,m);
end

%% Image and Reference Composite
x = [img, ref];

%% Squared Fourier Transform Magnitude 
xpad = zeros(m);
xpad(1:n, 1:2*n) = x;
f = fft2(xpad);
y_clean = abs(f).^2; 

%% Iterates
for photon_param_ind=1:length(photon_param_range)
    photon_param=photon_param_range(photon_param_ind);
    for trial=1:num_trials
        [photon_param,trial]
        %% Add noise
        n_photon = photon_param*m^2;
        y = sum(y_clean(:))*(n_photon)^(-1) * poissrnd( (n_photon/sum(y_clean(:)))*y_clean );
        %% Get autocorrelation
        rfull=(1/(m^2))*(F1c')*y*transpose((F2c'));
        r = rfull(1:n,1:n);
        %% Run Referenced Deconvolution Algorithm
        z=img_recov(r, n, ref, ref_type);
        %% Observed Error
        err_iter_vec(trial)=(norm(img(:)-z(:))/norm(img(:)))^2;
    end
    %% Average Observed Error
    err=sum(err_iter_vec)/num_trials;
    %% Expected Error
    c=n_photon/sum(y_clean(:));
    exp_err=(1/c)*sum(S(:).*y_clean(:))/(norm(img(:))^2);
    err_vec(photon_param_ind)=err;
    exp_err_vec(photon_param_ind)=exp_err;
    %% Save Recovered Image (for specified photons per pixel value)
    if photon_param_range(photon_param_ind)==photon_param_rec
        img_rec=real(z);
    end
end