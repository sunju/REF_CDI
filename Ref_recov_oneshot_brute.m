clear all
rng(1)
%% Generate problem size data
nimg = [150,150]; % Total x and y pixels in original image
nref = [150,150]; % Total x and y pixels in adjacent reference
% Note: nimg(1) and nref(1) need to be the same size in this implementation
len = [1001,1001]; % Total x and y pixels in FFT data (needs to be sufficiently large and ODD) 
cen = [5,5]; % Total number of x and y center pixels missing (needs to be ODD)
thresh=1e-99; % Threshold for smallest singular value to no be discarded
noise_model = 1; % 0 = no noise, 1 = unscaled Poisson, 2 = scaled Poisson (from Ponan)
solver = 1; % 0 = ifft implementation (assumes no missing center), 1 = least squares solver (allows for missing center)
% Input specimen
namestr = 'lena' ;
stanstr = 'jpg'      ;
X = mat2gray(imread([namestr,'.',stanstr])) ;
X_0 = rgb2gray(X);
img=imresize(X_0,nimg);
%img=ones(nimg);
%
tic
% Block referenc
block_ref=ones(nref);
ref=block_ref;
ref_type='b';
nimg=size(img);
nref=size(ref);
k=cen;
alpha=1;

x = [img, ref];
n = size(x);
xpad = zeros(len);
xpad(1:n(1), 1:n(2)) = x;
f = fft2(xpad);
y = abs(f).^2; 
% Add noise
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
% Get autocorrelation
% Remove center frequencies from data
four_supp1=1:len(1);
four_supp2=1:len(2);
if cen(1)>0
    four_supp1 = 1+(k(1)-1)/2+1:len(1) - (k(1)-1)/2;
end
if cen(2)>0
    four_supp2 = 1+(k(2)-1)/2+1:len(2) - (k(2)-1)/2;
end
yt = y(four_supp1, four_supp2);
%% Get autocorrelation
F1=dftmtx(len(1));
F1c=F1(:,[end-(n(1)-1)+1:end,1:n(1)]);
if cen(1)>0
    F1c=F1c(1+(cen(1)-1)/2+1:end-(cen(1)-1)/2,:);
end
[U1,D1,V1]=svd(F1c);
d1=diag(D1);
for t=1:length(d1)
    if d1(t)<thresh
        break;
    end
end
D1=D1(:,1:t);
V1=V1(:,1:t);
pinvF1c=V1*pinv(D1)*U1';
F2=dftmtx(len(2));
F2c=F2(:,[end-(n(2)-1)+1:end,1:n(2)]);
if cen(2)>0
    F2c=F2c(1+(cen(2)-1)/2+1:end-(cen(2)-1)/2,:);
end
[U2,D2,V2]=svd(F2c);
d2=diag(D2);
for t=1:length(d2)
    if d2(t)<thresh
        break;
    end
end
D2=D2(:,1:t);
V2=V2(:,1:t);
pinvF2c=V2*pinv(D2)*U2';
%% Get autocorrelation
rfull=pinvF1c*yt*transpose(pinvF2c);
r = real(rfull(1:nimg(1),1:nimg(2)));
%% Run algorithm
z=img_recov(r, nimg, alpha, ref, ref_type);
%% Observed error
err=norm(img(:)-z(:))/norm(img(:));
%%
toc
imshow(z)
err