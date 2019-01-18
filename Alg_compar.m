clear all
%% Generate problem input parameters
n = 64; % image and reference dimension
m = 1024; % detector dimension 
photon_param_range=[100,500,1000,2500,5000]; % range of photons per pixel values to be simulated
photon_param_rec=[1000]; % one chosen photon per pixel value for recovered image
ref = []; % specify reference if ref_type='a'; otherwise can leave empty
num_trials = 1; % number of trials (for each photon per pixel value) to average over
%% Generate problem input specimen
namestr = 'mimivirus' ;
stanstr = 'png'      ;
X = mat2gray(imread([namestr,'.',stanstr])) ;
X_0 = rgb2gray(X);
img=imresize(X_0,[n,n]);
%% Run algorithm
[block_err_vec,block_exp_err_vec,block_img_rec]=Ref_Deconv_fcn(n,m,img,'b',ref,num_trials,photon_param_range,photon_param_rec);
[slit_err_vec,slit_exp_err_vec,slit_img_rec]=Ref_Deconv_fcn(n,m,img,'s',ref,num_trials,photon_param_range,photon_param_rec);
[pinhole_err_vec,pinhole_exp_err_vec,pinhole_img_rec]=Ref_Deconv_fcn(n,m,img,'p',ref,num_trials,photon_param_range,photon_param_rec);
[HIO_block_err_vec,HIO_block_img_rec]=Ref_HIO_fcn(n,m,img,'b',ref,num_trials,20,photon_param_range,photon_param_rec);
[HIO_slit_err_vec,HIO_slit_img_rec]=Ref_HIO_fcn(n,m,img,'s',ref,num_trials,20,photon_param_range,photon_param_rec);
[HIO_pinhole_err_vec,HIO_pinhole_img_rec]=Ref_HIO_fcn(n,m,img,'p',ref,num_trials,20,photon_param_range,photon_param_rec);
[HIO_err_vec,HIO_img_rec]=HIO_fcn(n,m,img,num_trials,20,photon_param_range,photon_param_rec);
%%
save('Alg_compar.mat')
%%
figure;
semilogy(photon_param_range,block_err_vec,'red','LineWidth',1,'LineStyle','none','Marker','o')
hold on;
semilogy(photon_param_range,block_exp_err_vec,'red','LineWidth',1,'LineStyle','none','Marker','+')
hold on;
semilogy(photon_param_range,slit_err_vec,'blue','LineWidth',1,'LineStyle','none','Marker','o')
hold on;
semilogy(photon_param_range,slit_exp_err_vec,'blue','LineWidth',1,'LineStyle','none','Marker','+')
hold on;
semilogy(photon_param_range,pinhole_err_vec,'magenta','LineWidth',1,'LineStyle','none','Marker','o')
hold on;
semilogy(photon_param_range,pinhole_exp_err_vec,'magenta','LineWidth',1,'LineStyle','none','Marker','+')
hold on;
semilogy(photon_param_range,HIO_block_err_vec,'red','LineWidth',1,'LineStyle','none','Marker','*')
hold on
semilogy(photon_param_range,HIO_slit_err_vec,'blue','LineWidth',1,'LineStyle','none','Marker','*')
hold on
semilogy(photon_param_range,HIO_pinhole_err_vec,'magenta','LineWidth',1,'LineStyle','none','Marker','*')
hold on
semilogy(photon_param_range,HIO_err_vec,'black','LineWidth',1,'LineStyle','none','Marker','*')
hold off
legend('Ref. Deconv. (Block) - Exp.','Ref. Deconv. (Block) - Theory','Ref. Deconv. (Slit) - Exp.','Ref. Deconv. (Slit) - Theory','Ref. Deconv. (Pinhole) - Exp.','Ref. Deconv. (Pinhole) - Theory','HIO (Block)','HIO (Slit)','HIO (Pinhole)','HIO (No ref.)')
xlabel('Photons per pixel')
%ylim([0,0.007])