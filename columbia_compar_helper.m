clear all
%% Generate problem input parameters
n = 64; % image and reference dimension
m = 1024; % detector dimension 
photon_param_range=[1000]; % range of photons per pixel values to be simulated
photon_param_rec=[1000]; % one chosen photon per pixel value for recovered image
ref = []; % specify reference if ref_type='a'; otherwise can leave empty
num_trials = 1; % number of trials (for each photon per pixel value) to average over
%% Generate problem input specimen
namestr = 'columbia' ;
stanstr = 'jpg'      ;
X = mat2gray(imread([namestr,'.',stanstr])) ;
X_0 = rgb2gray(X);
img=imresize(X_0,[n,n]);
%% Run algorithm
[block_err_vec,block_exp_err_vec,block_img_rec]=Ref_Deconv_fcn(n,m,img,'b',ref,num_trials,photon_param_range,photon_param_rec);
[HIO_block_err_vec,HIO_block_img_rec]=Ref_HIO_fcn(n,m,img,'b',ref,num_trials,20,photon_param_range,photon_param_rec);
[HIO_err_vec,HIO_img_rec]=HIO_fcn(n,m,img,num_trials,20,photon_param_range,photon_param_rec);
%%
save('columbia_compar_helper.mat')
imwrite(block_img_rec,'block_columbia.png');
imwrite(HIO_block_img_rec,'HIO_block_columbia.png');
imwrite(HIO_img_rec,'HIO_columbia.png');