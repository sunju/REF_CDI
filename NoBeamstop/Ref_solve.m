clear all
rng(1)
%% Generate problem size data
nimg = [150,150]; % Total x and y pixels in original image
nref = [150,150]; % Total x and y pixels in adjacent reference
% Note: nimg(1) and nref(1) need to be the same size in this implementation
len = [1001,2001]; % Total x and y pixels in FFT data (needs to be sufficiently large and ODD) 
cen = [15,15];
ref_choice='b'; %Choose from 'b' (block), 'h' (holog), 'L' (L-shape), 'r' (random-Bernoulli)
noise_model = 2; % 0 = no noise, 1 = unscaled Poisson, 2 = scaled Poisson (from Ponan)
solver = 1; % 0 = no beamstop, 1 = beamstop
namestr = 'mimivirus' ;
stanstr = 'png'      ;
X = mat2gray(imread([namestr,'.',stanstr])) ;
X_0 = rgb2gray(X);
img=imresize(X_0,nimg);
%img=ones(nimg);
%
helper_mtrx=eye(nref(1)^2);

if ref_choice=='b'
    block_ref=ones(nref);
    ref=block_ref;
    ref_type='b';
end
if ref_choice=='h'
    holog_ref=zeros(nref); holog_ref(nref,nref)=1;
    ref=holog_ref;
    ref_type='h';
end
if ref_choice=='L'
    L_ref=zeros(nref); L_ref(:,nref(2))=ones(nref(2),1); L_ref(nref(1),:)=ones(1,nref(1));
    ref=L_ref;
    ref_type='L';
    helper_mtrx=sparse(ref2mtrx(ref));
end
if ref_choice=='r'
    r_ref=rand(nref)>0.5;
    r_ref(end,end)=1;
    ref=r_ref;
    ref_type='r';
    helper_mtrx=sparse(ref2mtrx(ref));
end
nimg=size(img);
nref=size(ref);

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
if noise_model == 3
    n_photon_order = 9;
    n_photon = 1.67 * 10^n_photon_order;
    nor_fac = norm(f(:),1);
    f = nor_fac * sqrt(  n_photon^-1 * poissrnd( n_photon/nor_fac^2 * abs(f).^2 ));
    y = f.^2;
end
%% Remove center frequencies from data
if solver == 1
    four_supp1=1:len(1);
    four_supp2=1:len(2);
    k=cen;
    if cen(1)>0
        four_supp1 = 1+(k(1)-1)/2+1:len(1) - (k(1)-1)/2;
    end
    if cen(2)>0
        four_supp2 = 1+(k(2)-1)/2+1:len(2) - (k(2)-1)/2;
    end
    yt = y(four_supp1, four_supp2);
    F1=dftmtx(len(1));
    F1c=F1(:,[end-(n(1)-1)+1:end,1:n(1)]);
    if cen(1)>0
        F1c=F1c(1+(cen(1)-1)/2+1:end-(cen(1)-1)/2,:);
    end
    F2=dftmtx(len(2));
    F2c=F2(:,[end-(n(2)-1)+1:end,1:n(2)]);
    if cen(2)>0
        F2c=F2c(1+(cen(2)-1)/2+1:end-(cen(2)-1)/2,:);
    end
    pinvF1c=pinv(F1c);
    pinvF2c=pinv(F2c);
end
%%
tic;
%% Get autocorrelation
if solver == 0
    r = ifft2(y);
    r = r([end-(n(1)-1)+1:end,1:n(1)],[end-(n(2)-1)+1:end,1:n(2)]);
    r = real(r(1:nimg(1),1:nimg(2)));
end
if solver == 1
    rfull=pinvF1c*yt*transpose(pinvF2c);
    r = real(rfull(1:nimg(1),1:nimg(2)));
end
%% Run algorithm
z=img_recov(r, nimg, ref, ref_type, helper_mtrx);
%%
toc;
%% Observed error
err=norm(img(:)-z(:))/norm(img(:));
%%
imshow(z)
err