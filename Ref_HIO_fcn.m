function [err_vec,img_rec]=Ref_HIO_fcn(n,m,img0,ref_type,ref,num_trials,max_iter,photon_param_range,photon_param_rec)

rng(1)
beta=0.5;

%% Reference choice

if ref_type=='b';
    ref=ones(n);
end

if ref_type=='s';
    ref=zeros(n); ref(:,end)=ones(n,1);
end

if ref_type=='p';
    ref=zeros(n); ref(n,n)=1;
end

%% Squared Fourier Magnitudes

X0=[img0,ref];
X0_os=zeros(m,m);
X0_os(m/2-n/2+1:m/2+n/2,m/2-2*n/2+1:m/2+2*n/2)=X0;
Y0_clean=sqrt(abs(fft2(X0_os)).^2);

%% Iterates
for photon_param_ind=1:length(photon_param_range)
    photon_param=photon_param_range(photon_param_ind);
    for trial=1:num_trials
        %% Add Noise
        y0_clean=Y0_clean.^2;
        n_photon = photon_param*m^2;
        y0 = sum(y0_clean(:))*(n_photon)^(-1) * poissrnd( (n_photon/sum(y0_clean(:)))*y0_clean );
        Y0=sqrt(y0);
        
        %% Initialize
        img_init=rand(n,n);
        X_init=[img_init,ref];
        X_os_init=zeros(m,m);
        X_os_init(m/2-n/2+1:m/2+n/2,m/2-2*n/2+1:m/2+2*n/2)=X_init;
        img_prev_init=rand(n,n);
        X_prev_init=[img_prev_init,ref];
        X_os_prev_init=zeros(m,m);
        X_os_prev_init(m/2-n/2+1:m/2+n/2,m/2-2*n/2+1:m/2+2*n/2)=X_prev_init;
        img=img_init;
        X_os=X_os_init;
        X_os_prev=X_os_prev_init;
        %% HIO Iterates
        for i=1:max_iter
            [photon_param,trial,i,norm(img0-img,'fro')/norm(img0,'fro')]
            Y=fft2(X_os);
            Y(1:m/2,:)=Y0(1:m/2,:).*(Y(1:m/2,:)./abs(Y(1:m/2,:)));
            Y(m/2+1:end,:)=Y0(m/2+1:end,:).*(Y(m/2+1:end,:)./abs(Y(m/2+1:end,:)));
            Y(m/2+1:m/2,1:m/2)=Y0(m/2+1:m/2,1:m/2).*(Y(m/2+1:m/2,1:m/2)./abs(Y(m/2+1:m/2,1:m/2)));
            Y(m/2+1:m/2,m/2+1:end)=Y0(m/2+1:m/2,m/2+1:end).*(Y(m/2+1:m/2,m/2+1:end)./abs(Y(m/2+1:m/2,m/2+1:end)));

            X_os=ifft2(Y);
            X_os(1:m/2-n/2,:) = X_os_prev(1:m/2-n/2,:) - beta*X_os(1:m/2-n/2,:);
            X_os(m/2+n/2+1:end,:) = X_os_prev(m/2+n/2+1:end,:) - beta*X_os(m/2+n/2+1:end,:);
            X_os(m/2-n/2+1:m/2+n/2,1:m/2-2*n/2) = X_os_prev(m/2-n/2+1:m/2+n/2,1:m/2-2*n/2) - beta*X_os(m/2-n/2+1:m/2+n/2,1:m/2-2*n/2);
            X_os(m/2-n/2+1:m/2+n/2,m/2+2*n/2+1:end) = X_os_prev(m/2-n/2+1:m/2+n/2,m/2+2*n/2+1:end) - beta*X_os(m/2-n/2+1:m/2+n/2,m/2+2*n/2+1:end);
            %% Enforce Known Reference
            X_os(m/2-n/2+1:m/2+n/2,m/2+1:m/2+2*n/2)=ref;
            %%
            X_os_prev=X_os;
            img=X_os(m/2-n/2+1:m/2+n/2,m/2-2*n/2+1:m/2+2*n/2);
            img=X_os(m/2-n/2+1:m/2+n/2,m/2-2*n/2+1:m/2);
        end
    err_iter_vec(trial)=(norm(img0(:)-img(:))/norm(img0(:)))^2;
    if photon_param_range(photon_param_ind)==photon_param_rec
        img_rec=real(img);
    end
    end
    err=sum(err_iter_vec)/num_trials;
    err_vec(photon_param_ind)=err;
end