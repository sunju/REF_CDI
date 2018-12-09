function z=img_recov(r, n, ref, ref_type);

%% Run algorithm
if ref_type=='b'
        k=tril(ones(n));
        z = k \ r / k';
end

if ref_type=='s'
        k=tril(ones(n));
        z = k \ r / eye(n)';
end

if ref_type=='p'
            z=r;            
end
if ref_type=='a'
    L=sparse(ref2mtrx(ref));
    b=r; b=b(:);
    z=L\b;
end
z=reshape(z,n,n);