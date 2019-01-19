% Reproduces the spectrum in Figure 9 of the paper
% "Holographic Phase Retrieval and Optimal Reference Design",
% by David Barmherzig, Ju Sun, Emmanuel Candes, T.J. Lane, and Po-Nan Li.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clear all
rng(1)
n=64;
m=1024;
X = imread('mimivirus.png');
X = mat2gray(X);
X_0 = rgb2gray(X);
img=imresize(X_0,[n,n]);
imgpad=zeros(m,m); imgpad(1:n,1:n)=img;
y=abs(fft2(imgpad)).^2;
%%
figure;
stem3([-512:511],[-512:511],fftshift(y))