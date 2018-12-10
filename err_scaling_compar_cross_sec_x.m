clear all
n=64;
m=1024;
bndry=511;
for k1=-bndry:bndry+1
    k2=0;
    lower_bnd(k1+bndry+1)=1/(m^4);
    scaler_b(k1+bndry+1)=( (1/(m^2))+2*((n-1)/(m^2))*(1-cos(2*pi*k1/m)) )*( (1/(m^2))+2*((n-1)/(m^2))*(1-cos(2*pi*k2/m)));
    scaler_p(k1+bndry+1)=((n^2)/(m^4));
    scaler_s(k1+bndry+1)=(n/(m^2))*( (1/(m^2))+2*((n-1)/(m^2))*(1-cos(2*pi*k2/m)));
end
figure;
plot([-bndry:bndry+1],lower_bnd,'red','LineWidth',2)
hold on;
plot([-bndry:bndry+1],scaler_b,'blue','LineWidth',2)
hold on;
plot([-bndry:bndry+1],scaler_p,'green','LineWidth',2)
hold on;
plot([-bndry:bndry+1],scaler_s,'magenta','LineWidth',2)
hold off
legend({'Uniform lower bound','Block ref.', 'Pinhole ref.', 'Slit ref.'},'FontSize',14)
xlabel('x-frequency','FontSize',18)
xlim([-512,512])