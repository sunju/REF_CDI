function x=fft_handle_solver2d(b, n, len, four_supp1, four_supp2)

tol = 1e-8; 
maxit = 2000;

four_supp_row1 = four_supp1 - 1;
four_supp_col1 = [len(1) - (n(1) - 2) : len(1), 1 : n(1)] - 1;
four_supp_row2 = four_supp2 - 1;
four_supp_col2 = [len(2) - (n(2) - 2) : len(2), 1 : n(2)] - 1;

x=lsqr(@fft_handle,b,tol,maxit);

    function y = fft_handle(x,transp_flag)
        if strcmp(transp_flag,'transp')      % y = A'*x
            y=exp(2*pi*i*four_supp_col1'*four_supp_row1/len(1))*reshape(x,length(four_supp_row1),length(four_supp_row2))*exp(2*pi*i*four_supp_row2'*four_supp_col2/len(2));
            y=vec(y);
        elseif strcmp(transp_flag,'notransp') % y = A*x
            y=exp(-2*pi*i*four_supp_row1'*four_supp_col1/len(1))*reshape(x,length(four_supp_col1),length(four_supp_col2))*exp(-2*pi*i*four_supp_col2'*four_supp_row2/len(2));
            y=vec(y);
        end
    end

end

%kron(transpose(F2st),F1st)