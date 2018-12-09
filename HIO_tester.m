%% Generate problem input parameters
n = 64; % image and reference dimension
m = 1024; % detector dimension 
photon_param_range=[100,500,1000,2500,5000]; % range of photons per pixel values to be simulated
photon_param_rec=[1000]; % one chosen photon per pixel value for recovered image
num_trials=1; % number of trials (for each photon per pixel value) to average over
max_iter=15;
%% Generate problem input specimen
namestr = 'mimivirus' ;
stanstr = 'png'      ;
X = mat2gray(imread([namestr,'.',stanstr])) ;
X_0 = rgb2gray(X);
img0=imresize(X_0,[n,n]);
%% Run algorithm
[err_vec,img_rec]=HIO_fcn(n,m,img0,num_trials,max_iter,photon_param_range,photon_param_rec);
err_vec
imshow(img_rec)
save('Ref_HIO_tester.mat')