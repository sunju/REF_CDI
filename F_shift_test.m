clear all
n = 10;
img = rand(n,1);
x1 = [img; zeros(n,1)];
x2 = [zeros(n,1); img];
y1 = abs(fft(x1)).^2;
y2 = abs(fft(x2)).^2;
norm(y1-y2)