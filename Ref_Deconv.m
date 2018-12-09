clear all
rng(1)
%% Generate problem size data
n = 64; % Image and reference dimension
m = 1024; % Detector dimension 
photon_param_range=[100,500,1000,2500,5000];
ref_type = 's'; %'b'=block, 's'=slit, 'p'=pinhole
num_trials=1;
% Input specimen
namestr = 'mimivirus' ;
stanstr = 'png'      ;
X = mat2gray(imread([namestr,'.',stanstr])) ;
X_0 = rgb2gray(X);
img=imresize(X_0,[n,n]);

block_ref=ones(n);
slit_ref=zeros(n); slit_ref(:,end)=ones(n,1);
pinhole_ref=zeros(n); pinhole_ref(n,n)=1;

F1=dftmtx(m);
F1=F1(:,[end-(n-1)+1:end,1:n]);
pinvF1=pinv(F1);
F2=dftmtx(m);
F2=F2(:,[end-(2*n-1)+1:end,1:2*n]);
pinvF2=pinv(F2);

Dn=diffmtrx(n);

S_flat = (n/(m^2))*ones(1,m);
S_dec = 1/(m^2) + 2*((n-1)/(m^2))*(1-cos(2*pi*[0:m-1]/m));

if ref_type=='b';
    ref=block_ref;
    MRinv1=Dn;
    MRinv2=Dn;
    S=S_dec'*S_dec;
end

if ref_type=='p';
    ref=pinhole_ref;
    MRinv1=eye(n);
    MRinv2=eye(n);
    S=S_flat'*S_flat;
end

if ref_type=='s';
    ref=slit_ref;
    MRinv1=eye(n);
    MRinv2=Dn;
    S=S_dec'*S_flat;
end

x = [img, ref];

xpad = zeros(m);
xpad(1:n, 1:2*n) = x;
f = fft2(xpad);
y_clean = abs(f).^2; 
% Add noise
% Get autocorrelation
for photon_param_ind=1:length(photon_param_range)
    photon_param=photon_param_range(photon_param_ind);
    for trial=1:num_trials
        [photon_param,trial]
        %% Add noise
        n_photon = photon_param*m^2;
        y = sum(y_clean(:))*(n_photon)^(-1) * poissrnd( (n_photon/sum(y_clean(:)))*y_clean );
        %% Get autocorrelation
        rfull=pinvF1*y*transpose(pinvF2);
        r = rfull(1:n,1:n);
        %% Run algorithm
        z=img_recov(r, n, ref, ref_type);
        %% Observed error
        err_iter_vec(trial)=norm(img(:)-z(:))/norm(img(:));
    end
    err=sum(err_iter_vec)/num_trials;
    %% Expected error
    c=n_photon/sum(y_clean(:));
    exp_err1=sqrt((1/c)*sum(T(:).*y_clean(:)))/norm(img(:));
    err_vec(photon_param_ind)=err;
    exp_err_vec(photon_param_ind)=exp_err;
end
err_vec
exp_err
save('Ref_Deconv.mat')