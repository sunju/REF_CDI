% solve the beam stopping LS for 1D case
n = 100;
x = randn(n, 1);     % true signal
x_corr = conv(flip(x), x);

over_sample = 4;
obs_length = round(n*over_sample);
BS_fact = 0.04;
BS_num = round(BS_fact*obs_length);
y = abs(fft(x, n*over_sample)).^2;
% this reproduces x_corr with zeros padded at two ends
fftshift(ifft(y));

 % for simplicity, we assume y is always of even-length, as well as BS_num
BS_y = fftshift(y);
BS_y((length(y)/2-BS_num/2+1):(length(y)/2+BS_num/2)) = [];
length(BS_y)
