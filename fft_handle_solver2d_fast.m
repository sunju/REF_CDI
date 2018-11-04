function x=fft_handle_solver2d_fast(tol, maxit, b, n, len, four_supp1, four_supp2, cen)

x=lsqr(@fft_handle,b,tol,maxit);

    function y = fft_handle(x,transp_flag)
        if strcmp(transp_flag,'transp')      % y = A'*x
            x = reshape(x,len(1)-cen(1),len(2)-cen(2));
            x_cen_adj=zeros(len(1),len(2));
            x_cen_adj(1+(cen(1)-1)/2+1:1+(cen(1)-1)/2+(len(1)-cen(1)),1+(cen(2)-1)/2+1:1+(cen(2)-1)/2+(len(2)-cen(2)))=x;
            y=len(1)*len(2)*ifft2(x_cen_adj);
            y=circshift(y,[n(1)-1,n(2)-1]);
            y=y(1:2*n(1)-1,1:2*n(2)-1);
            y=real(vec(y));
        elseif strcmp(transp_flag,'notransp') % y = A*x
            xps=zeros(len(1),len(2)); xps(1:2*n(1)-1,1:2*n(2)-1)=reshape(x,2*n(1)-1,2*n(2)-1);
            xps=circshift(xps,[-(n(1)-1),-(n(2)-1)]);
            y=fft2(xps);
            if (cen(1)>0)&(cen(2)>0)
                y=y(1+(cen(1)-1)/2+1:end-(cen(1)-1)/2,1+(cen(2)-1)/2+1:end-(cen(2)-1)/2);
            end
            y=real(vec(y));
        end
    end
end