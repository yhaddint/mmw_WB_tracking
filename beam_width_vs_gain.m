clear;clc;
N = 2e3;
xdata = linspace(-pi/2,pi/2,N);
mu = 0;
b = 2.5/180*pi;

lappdf = 1/(2*b)*exp(-abs(xdata-mu)/b);
figure
plot(xdata,lappdf);hold on
grid on

beamwidth_num = 256;
beamwidth_range = linspace(pi/beamwidth_num,pi,beamwidth_num-1);
for bb = 1:beamwidth_num/4
    beam_power = double(abs(xdata) - beamwidth_range(bb)<0);
    beam_pattern = beam_power./norm(beam_power);
    gain(bb) = sum(lappdf.*beam_pattern);
end

figure
plot(beamwidth_range(1:64)/pi*180,10*log10(gain))
grid on
xlabel('beamwidth (deg)')
ylabel('Gain (dB)')
