clear all
n=2;
m=4;
x=rand(n);
%% Three ways to get the oversampled squared FT
xpad=zeros(m,m); xpad(1:n,1:n)=x;
y1=abs(fft2(xpad)).^2
%%
y2=abs(fft2(x,m,m)).^2
F=dftmtx(m);
y3=abs(F(:,1:n)*x*transpose(F(:,1:n))).^2
%% Autocorrelation
ax=conv2(x,flip(flip(x,1),2))%%
ax_pad=zeros(m,m); ax_pad(1:2*n-1,1:2*n-1)=ax;
% Shift autocorrelation to get center value at (0,0) index
ax_pad_s=circshift(circshift(ax_pad,-(n-1),1),-(n-1),2); % note this is needed!
y4=F*ax_pad_s*transpose(F) %% equals y !
ax_s=circshift(circshift(ax,-(n-1),1),-(n-1),2);
y5=F(:,[1:n,end-(n-1)+1:end])*ax_s*transpose(F(:,[1:n,end-(n-1)+1:end]))
norm(y5-y1)
