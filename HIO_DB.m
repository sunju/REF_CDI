clear all
rng(1)
n=150; % dimensions of the image (and reference if used)
L=1024; % dimensions of the data
c=[50,50]; % dimensions of the missing center square
beta=0.5;
max_iter=200;
err_thresh=0.05;
ref_flag=1; % set to 0 if no reference, and to 1 if reference added to image
ref_constr_flag=1; % (applicable when ref_flag==1), set to 1 to impose reference pixels during HIO algorithm
%%
namestr = 'mimivirus' ;
stanstr = 'png'      ;
img0       = mat2gray(imread([namestr,'.',stanstr])) ;
img0 = rgb2gray(img0);
img0=imresize(img0,[n,n]);
%%
if ref_flag==1
    %img0=rand(n,n);
    %ref=zeros(n,n); ref(n,1:n)=ones(1,n); ref(1:n,n)=ones(n,1); %L-shape
    ref=ones(n,n); %Block-ref
    %ref=zeros(n,n); ref(n,n)=1; %Point ref
end
%%
X0=img0;
if ref_flag==1
    X0=[img0,ref];
end
n1=size(X0,1); n2=size(X0,2);
X0_os=zeros(L,L);
X0_os(L/2-n1/2+1:L/2+n1/2,L/2-n2/2+1:L/2+n2/2)=X0;
Y0=sqrt(poissrnd(abs(fft2(X0_os)).^2));
%%
img_init=rand(n,n);
X_init=img_init;
if ref_flag==1
    X_init=[img_init,ref];
end
X_os_init=zeros(L,L);
X_os_init(L/2-n1/2+1:L/2+n1/2,L/2-n2/2+1:L/2+n2/2)=X_init;
%%
img_prev_init=rand(n,n);
X_prev_init=img_prev_init;
if ref_flag==1
    X_prev_init=[img_prev_init,ref];
end
X_os_prev_init=zeros(L,L);
X_os_prev_init(L/2-n1/2+1:L/2+n1/2,L/2-n2/2+1:L/2+n2/2)=X_prev_init;

img=img_init;
X_os=X_os_init;
X_os_prev=X_os_prev_init;

for i=1:max_iter
    [i, norm(img0-img,'fro')/norm(img0,'fro')]
    imshow(img)
    Y=fft2(X_os);
    Y(1:L/2-c/2,:)=Y0(1:L/2-c/2,:).*(Y(1:L/2-c/2,:)./abs(Y(1:L/2-c/2,:)));
    Y(L/2+c/2+1:end,:)=Y0(L/2+c/2+1:end,:).*(Y(L/2+c/2+1:end,:)./abs(Y(L/2+c/2+1:end,:)));
    Y(L/2-c/2+1:L/2+c/2,1:L/2-c/2)=Y0(L/2-c/2+1:L/2+c/2,1:L/2-c/2).*(Y(L/2-c/2+1:L/2+c/2,1:L/2-c/2)./abs(Y(L/2-c/2+1:L/2+c/2,1:L/2-c/2)));
    Y(L/2-c/2+1:L/2+c/2,L/2+c/2+1:end)=Y0(L/2-c/2+1:L/2+c/2,L/2+c/2+1:end).*(Y(L/2-c/2+1:L/2+c/2,L/2+c/2+1:end)./abs(Y(L/2-c/2+1:L/2+c/2,L/2+c/2+1:end)));
    Y(L/2-c/2+1:L/2+c/2,L/2-c/2+1:L/2+c/2)=Y0(L/2-c/2+1:L/2+c/2,L/2-c/2+1:L/2+c/2).*(Y(L/2-c/2+1:L/2+c/2,L/2-c/2+1:L/2+c/2)./abs(Y(L/2-c/2+1:L/2+c/2,L/2-c/2+1:L/2+c/2)));
    %Y=Y0.*(Y./abs(Y));
    X_os=ifft2(Y);
    X_os(1:L/2-n1/2,:) = X_os_prev(1:L/2-n1/2,:) - beta*X_os(1:L/2-n1/2,:);
    X_os(L/2+n1/2+1:end,:) = X_os_prev(L/2+n1/2+1:end,:) - beta*X_os(L/2+n1/2+1:end,:);
    X_os(L/2-n1/2+1:L/2+n1/2,1:L/2-n2/2) = X_os_prev(L/2-n1/2+1:L/2+n1/2,1:L/2-n2/2) - beta*X_os(L/2-n1/2+1:L/2+n1/2,1:L/2-n2/2);
    X_os(L/2-n1/2+1:L/2+n1/2,L/2+n2/2+1:end) = X_os_prev(L/2-n1/2+1:L/2+n1/2,L/2+n2/2+1:end) - beta*X_os(L/2-n1/2+1:L/2+n1/2,L/2+n2/2+1:end);
    if (ref_flag==1)&(ref_constr_flag==1)
        X_os(L/2-n1/2+1:L/2+n1/2,L/2+1:L/2+n2/2)=ref;
    end
    X_os_prev=X_os;
    img=X_os(L/2-n1/2+1:L/2+n1/2,L/2-n2/2+1:L/2+n2/2);
    if ref_flag==1
        img=X_os(L/2-n1/2+1:L/2+n1/2,L/2-n2/2+1:L/2);
    end
    if norm(img-img0,'fro')/norm(img0,'fro')<err_thresh
        break;
    end
end
i, norm(img0-img,'fro')/norm(img0,'fro')