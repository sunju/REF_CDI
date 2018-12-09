clear all
n=64;
m=1024;
bndry=500;
for k1=-bndry:bndry
    for k2=-bndry:bndry
        lower_bnd(k1+bndry+1,k2+bndry+1)=1/(m^4);
        scaler_b(k1+bndry+1,k2+bndry+1)=( (1/(m^2))+2*((n-1)/(m^2))*(1-cos(2*pi*k1/m)) )*( (1/(m^2))+2*((n-1)/(m^2))*(1-cos(2*pi*k2/m)));
        scaler_p(k1+bndry+1,k2+bndry+1)=((n^2)/(m^4));
        scaler_s(k1+bndry+1,k2+bndry+1)=(n/(m^2))*( (1/(m^2))+2*((n-1)/(m^2))*(1-cos(2*pi*k2/m)));
    end
end
figure;
stem3([-bndry:bndry],[-bndry:bndry],lower_bnd,'red')
hold on;
stem3([-bndry:bndry],[-bndry:bndry],scaler_b,'blue')
hold on;
stem3([-bndry:bndry],[-bndry:bndry],scaler_p,'green')
hold on;
stem3([-bndry:bndry],[-bndry:bndry],scaler_s,'magenta')
hold off
%direction = [1 1 0];
%rotate(hSurface,direction,25)
legend('Uniform lower bound','Block ref.', 'Pinhole ref.', 'Slit ref.')
xlabel('x-frequency')
ylabel('y-frequency')
%ylim([0,1e-5])