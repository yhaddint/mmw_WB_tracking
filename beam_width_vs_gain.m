clear;clc;
N = 5e3;
xdata = linspace(-pi/2,pi/2,N);
mu = 0;
b = 4/180*pi;
% lappdf = 1/(2*b)*exp(-abs(xdata-mu)/b);
normalpdf = 1/sqrt(2*pi*b^2)*exp(-(xdata-mu).^2/(2*b^2));

% figure
% plot(xdata/pi*180,lappdf);
% xlabel('Angle Spread [deg]')
% ylabel('Probability')
% xlim([-30,30])
% hold on
% grid on
% %% sweep beamwidth
% steer_error = 5/180*pi;
% beamwidth_num = 2048;
% beamwidth_range = linspace(pi/beamwidth_num,pi,beamwidth_num-1);
% for bb = 1:beamwidth_num/8
%     Beam_BW = beamwidth_range(bb);
%     for nn = 1:N
%         inbeam = abs(xdata(nn)-steer_error)<=Beam_BW;
%         beam_env(nn) = double(inbeam);
%     end
%     beam_pow = pi/beamwidth_range(bb);
%     beam_pattern = beam_env * sqrt(beam_pow);
%     gain(bb) = sum(lappdf.*beam_pattern)*(xdata(2)-xdata(1));
% end
% 
% figure
% semilogx((beamwidth_range(1:256)/pi*180),20*log10(gain))
% grid on
% xlabel('Beamwidth [deg]')
% ylabel('Gain (dB)')
% xlim([1,16])
% ylim([0,25])

%% sweep estimation error
error_num = 1e3;
steer_error_range = linspace(0,20,error_num)/180*pi;


for ee = 1:error_num
    Beam_BW = 5/180*pi;
    steer_error = steer_error_range(ee);
    for nn = 1:N
        inbeam = abs(xdata(nn)-steer_error)<=Beam_BW/2;
        beam_env(nn) = double(inbeam);
    end
    beam_pow = pi/Beam_BW;
    beam_pattern = beam_env * sqrt(beam_pow);
    gain(ee) = sum(normalpdf.*beam_pattern)*(xdata(2)-xdata(1));
end

figure
plot(steer_error_range/pi*180,20*log10(gain))
grid on
xlabel('Steering Error [deg]')
ylabel('Gain (dB)')
xlim([1,16])
ylim([0,25])
