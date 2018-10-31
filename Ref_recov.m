clear all
rng(1)
%% Generate problem size data
nimg = [150,150]; % Total x and y pixels in original image
nref = [150, 150]; % Total x and y pixels in adjacent reference
% Note: nimg(1) and nref(1) need to be the same size in this implementation
len = [1024, 1024]; % Total x and y pixels in FFT data (needs to be sufficiently large) 
centersize = [1, 1]; % Total number of x and y center pixels missing
noise_model = 1; % 0 = no noise, 1 = unscaled Poisson, 2 = scaled Poisson (from Ponan)
solver = 1; % 0 = ifft implementation (assumes no missing center), 1 = least squares solver (allows for missing center)
%% Input specimen
namestr = 'mimivirus' ;
stanstr = 'png'      ;
X = mat2gray(imread([namestr,'.',stanstr])) ;
X_0 = rgb2gray(X);
img=imresize(X_0,nimg);
alpha_min=-5;
alpha_steps=25;
alpha_max=5;
num_trials=1;
%%
tic
%% Block reference
block_ref=ones(nref);
[block_alpha_vec,block_err_vec]=ref_solver(solver,noise_model,img,block_ref,'b',len,centersize,alpha_min,alpha_steps,alpha_max,num_trials);
%%% L-shaped reference
%L_ref=zeros(nref); L_ref(nref(1),:)=ones(1,nref(2)); L_ref(:,nref(2))=ones(nref(1),1);
%[L_alpha_vec,L_err_vec]=scale_trials_observed_nonunif(solver,noise_model,img,L_ref,'L',len,centersize,alpha_min,alpha_steps,alpha_max,num_trials);
%% Holography reference
holog_ref=zeros(nref); holog_ref(nref(1),nref(2))=1;
[holog_alpha_vec,holog_err_vec]=ref_solver(solver,noise_model,img,holog_ref,'h',len,centersize,alpha_min,alpha_steps,alpha_max,num_trials);
%%
toc
%% Plot
figure
loglog(block_alpha_vec,block_err_vec,'-xb')
hold on;
%loglog(L_alpha_vec,L_err_vec,'-xr')
%hold on;
loglog(holog_alpha_vec,holog_err_vec,'-xg')
%legend('Observed err. - Square ref.','Observed err. - L-shape ref.','Observed err. - Delta point ref.')
legend('Observed err. - Square ref.','Observed err. - Delta point ref.')
%title('Observed and expected relative error vs. reference scaling (100x100 mimvirus, with 10 trials per data point)')
xlabel('Reference-to-specimen transmittance ratio')
ylabel('Relative error')

save('Ref_recov.dat')