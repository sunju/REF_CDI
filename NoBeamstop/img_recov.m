function z=img_recov(r, nimg, ref, ref_type, helper_mtrx);

%% Run algorithm
if ref_type=='b'
        k=tril(ones(nimg(1),nimg(2)));
        z = k \ r / k';
end
if ref_type=='h'
            z=r(:);            
end
if (ref_type~='b')&(ref_type~='h')
    L=helper_mtrx;
    b=r; b=b(:);
    z=L\b;
end
z=reshape(z,nimg(1),nimg(2));
%z=real(reshape(z,n,n));