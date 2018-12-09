clear all
rng(1)
n=64;
m=1024;
X = imread('mimivirus.png');
X = mat2gray(X);
X_0 = rgb2gray(X);
img=imresize(X_0,[n,n]);
n=size(img,1); % assuming image is square
bref=ones(n,n);
pref=zeros(n,n); pref(end,end)=1;
sref=zeros(n,n); sref(:,end)=ones(n,1);
xb=[img,bref];
xbpad=zeros(m,m); xbpad(1:n,1:2*n)=xb;
yb=abs(fft2(xbpad)).^2;
xp=[img,pref];
xppad=zeros(m,m); xppad(1:n,1:2*n)=xp;
yp=abs(fft2(xppad)).^2;
xs=[img,sref];
xspad=zeros(m,m); xspad(1:n,1:2*n)=xs;
ys=abs(fft2(xspad)).^2;
%%
figure;
subplot(3,1,1)
stem3([-512:511],[-512:511],fftshift(yb),'red')
title('Squared Four. magnitudes for mimivirus with block ref.')
subplot(3,1,2)
stem3([-512:511],[-512:511],fftshift(ys),'blue')
title('Squared Four. magnitudes for mimivirus with slit ref.')
subplot(3,1,3)
stem3([-512:511],[-512:511],fftshift(yp),'magenta')
title('Squared Four. magnitudes for mimivirus with pinhole ref.')