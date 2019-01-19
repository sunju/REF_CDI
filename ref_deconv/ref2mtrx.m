% Function which maps an arbitrary reference ref to the corresponding
% matrix M_R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mtrx=ref2mtrx(ref)

[n1,n2]=size(ref); % assuming n1=n2
n=n1;

mtrx=zeros(n^2,n^2);
for i=1:n
    for j=1:i
        t=i-j+1;
        col=ref(:,n+1-t);
        blk=zeros(n,n);
        for k=1:n
            blk(k,1:k)=col(n+1-k:n);
        end
        mtrx(n*(i-1)+1:n*(i-1)+n,n*(j-1)+1:n*(j-1)+n)=blk;
    end
end