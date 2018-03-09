clear;clc;
N = 32;
angle_range = linspace(-pi,pi,2000);
for aa=1:length(angle_range)
    array_power(aa) = sum(exp(1j*pi*(0:N-1).'*sin(angle_range(aa))))/sqrt(N);
end

figure;plot(angle_range/pi*180,20*log10(abs(array_power)))
grid on
ylim([0,30])
%%
sum(abs(array_power).^2)*(angle_range(2)-angle_range(1))/(2*pi)