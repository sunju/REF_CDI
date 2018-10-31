function z=img_recov(r, nimg, alpha, ref, ref_type);

%% Run algorithm
if ref_type=='b'
        k=tril(sqrt(alpha)*ones(nimg(1),nimg(2)));
        z = k \ r / k';
end
if ref_type=='h'
            z=(1/alpha)*r(:);            
end
if (ref_type~='b')&(ref_type~='h')
    L=sparse(ref2mtrx(ref));
    b=r; b=b(:);
    z=L\b;
end
z=reshape(z,nimg(1),nimg(2));
%z=real(reshape(z,n,n));