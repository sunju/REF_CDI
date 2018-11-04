clear all
rng(1)
%% Generate problem size data
nimg = [150,150]; % Total x and y pixels in original image
nref = [150, 150]; % Total x and y pixels in adjacent reference
% Note: nimg(1) and nref(1) need to be the same size in this implementation
len = [1001, 1001]; % Total x and y pixels in FFT data (needs to be sufficiently large and ODD) 
centersize = [5,5]; % Total number of x and y center pixels missing (needs to be ODD)
noise_model = 1; % 0 = no noise, 1 = unscaled Poisson, 2 = scaled Poisson (from Ponan)
solver = 1; % 0 = ifft implementation (assumes no missing center), 1 = least squares solver (allows for missing center)
tol = 1e-8; 
maxit = 2000;
%% Input specimen
namestr = 'lena' ;
stanstr = 'jpg'      ;
X = mat2gray(imread([namestr,'.',stanstr])) ;
X_0 = rgb2gray(X);
img=imresize(X_0,nimg);
alpha_min=0; %-5;
alpha_steps= 1; %25;
alpha_max= 0; %5;
num_trials=1;
%%
tic
%% Block reference
block_ref=ones(nref);
[z,err]=ref_solver_oneshot_lsq(tol, maxit, solver,noise_model,img,block_ref,'b',len,centersize,alpha_min,alpha_steps,alpha_max,num_trials);
%% L-shaped reference
%L_ref=zeros(nref); L_ref(nref(1),:)=ones(1,nref(2)); L_ref(:,nref(2))=ones(nref(1),1);
%[z,err]=ref_solver_oneshot_lsq(tol, maxit, solver,noise_model,img,L_ref,'L',len,centersize,alpha_min,alpha_steps,alpha_max,num_trials);
%% Holography reference
%holog_ref=zeros(nref); holog_ref(nref(1),nref(2))=1;
%[z,err]=ref_solver_oneshot_lsq(tol, maxit, solver,noise_model,img,L_ref,'L',len,centersize,alpha_min,alpha_steps,alpha_max,num_trials);
%%
toc
imshow(z)
err