function Dn=diffmtrx(n);

Dn=eye(n);
for i=2:n
    Dn(i,i-1)=-1;
end