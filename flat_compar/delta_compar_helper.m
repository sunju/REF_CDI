% Reproduces the images shown in Figure 10 of the paper
% "Holographic Phase Retrieval and Optimal Reference Design",
% by David Barmherzig, Ju Sun, Emmanuel Candes, T.J. Lane, and Po-Nan Li.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% Generate problem input parameters
n = 64; % image and reference dimension
m = 1024; % detector dimension 
photon_param_range=[1000]; % range of photons per pixel values to be simulated
photon_param_rec=[1000]; % one chosen photon per pixel value for recovered image
ref = []; % specify reference if ref_type='a'; otherwise can leave empty
num_trials = 1; % number of trials (for each photon per pixel value) to average over
%% Generate problem input specimen
delta_rad=2;
img=zeros(n);
img(n/2-delta_rad:n/2+delta_rad,n/2-delta_rad:n/2+delta_rad)=ones(2*delta_rad+1,2*delta_rad+1);
%% Run algorithm
[block_err_vec,block_exp_err_vec,block_img_rec]=Ref_Deconv_fcn(n,m,img,'b',ref,num_trials,photon_param_range,photon_param_rec);
[slit_err_vec,slit_exp_err_vec,slit_img_rec]=Ref_Deconv_fcn(n,m,img,'s',ref,num_trials,photon_param_range,photon_param_rec);
[pinhole_err_vec,pinhole_exp_err_vec,pinhole_img_rec]=Ref_Deconv_fcn(n,m,img,'p',ref,num_trials,photon_param_range,photon_param_rec);
[HIO_block_err_vec,HIO_block_img_rec]=Ref_HIO_fcn(n,m,img,'b',ref,num_trials,20,photon_param_range,photon_param_rec);
[HIO_slit_err_vec,HIO_slit_img_rec]=Ref_HIO_fcn(n,m,img,'s',ref,num_trials,20,photon_param_range,photon_param_rec);
[HIO_pinhole_err_vec,HIO_pinhole_img_rec]=Ref_HIO_fcn(n,m,img,'p',ref,num_trials,20,photon_param_range,photon_param_rec);
[HIO_err_vec,HIO_img_rec]=HIO_fcn(n,m,img,num_trials,20,photon_param_range,photon_param_rec);
%%
save('delta_compar_helper.mat')
imwrite(block_img_rec,'block_delta.png');
imwrite(slit_img_rec,'slit_delta.png');
imwrite(pinhole_img_rec,'pinhole_delta.png');
imwrite(HIO_block_img_rec,'HIO_block_delta.png');
imwrite(HIO_slit_img_rec,'HIO_slit_delta.png');
imwrite(HIO_pinhole_img_rec,'HIO_pinhole_delta.png');
imwrite(HIO_img_rec,'HIO_delta.png');