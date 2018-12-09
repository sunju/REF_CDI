clear all
n=64;
m=1024;
bndry=500;
for k1=-bndry:bndry
    k2=0;
    lower_bnd(k1+bndry+1)=1/(m^4);
    scaler_b(k1+bndry+1)=( (1/(m^2))+2*((n-1)/(m^2))*(1-cos(2*pi*k1/m)) )*( (1/(m^2))+2*((n-1)/(m^2))*(1-cos(2*pi*k2/m)));
    scaler_p(k1+bndry+1)=((n^2)/(m^4));
    scaler_s(k1+bndry+1)=(n/(m^2))*( (1/(m^2))+2*((n-1)/(m^2))*(1-cos(2*pi*k2/m)));
end
figure;
plot([-bndry:bndry],lower_bnd,'red')
hold on;
plot([-bndry:bndry],scaler_b,'blue')
hold on;
plot([-bndry:bndry],scaler_p,'green')
hold on;
plot([-bndry:bndry],scaler_s,'magenta')
hold off
%direction = [1 1 0];
%rotate(hSurface,direction,25)
legend('Uniform lower bound','Block ref.', 'Pinhole ref.', 'Slit ref.')
xlabel('x-frequency')
%ylabel('y-frequency')
%ylim([0,1e-5])