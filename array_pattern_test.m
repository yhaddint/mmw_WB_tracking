clear;clc;
Nr = 32;
% %% Codebook - single beam
% beam_width = 60;
% angle_range = 180;
% pointing_range = -50:1:50;
% desired_pattern = zeros(angle_range*2,length(pointing_range));
% for ii=1:length(pointing_range)
%     left_boundary = pointing_range(ii) - beam_width/2;
%     right_boundary = pointing_range(ii) + beam_width/2 - 1;
%     index_desired = (left_boundary:right_boundary) + angle_range + 1;
%     desired_pattern(index_desired,ii) = ones(beam_width,1);
% end
% for kk = 1:angle_range*2
%     FF(:,kk) = exp(1j*pi*(0:Nr-1).'*sin((kk -1)/180*pi));
% end
% 
% for ii=1:101
%      w_data_raw = pinv(FF')*desired_pattern(:,ii);
%      w_data(:,ii) = w_data_raw./norm(w_data_raw);
% end
% %% Pattern Evaluation
% angle_test = -45;
% figure
% phi = (-180:179)/180*pi;
% % plot(phi/pi*180,20*log10(abs(FF'*w_data(:,angle_test+51))));
% polarplot(phi,(abs(FF'*w_data(:,angle_test+51))));
% hold on;
% % naive_BF_vec = exp(1j*pi*(0:Nr-1).'*sin(angle_test/180*pi))/sqrt(Nr);
% % polarplot(phi,20*log10(abs(FF'*naive_BF_vec)));
% % legend('codebook','naive')
% hold on
% grid on
% % ylim([-10,15])
%% Codebook - Four beam
beam_width = 4;
angle_range = 180;
pointing_range = -50:1:50;
desired_pattern = zeros(angle_range*2,length(pointing_range));
for ii=1:length(pointing_range)
    left_boundary = pointing_range(ii) - beam_width/2;
    right_boundary = pointing_range(ii) + beam_width/2 - 1;
    index_desired = (left_boundary:right_boundary) + angle_range + 1;
    desired_pattern(index_desired,ii) = ones(beam_width,1);
end
for kk = 1:angle_range*2
    FF(:,kk) = exp(1j*pi*(0:Nr-1).'*sin((kk -1)/180*pi));
end

angle_test = [-50, -30, -10, 50];

for ii=1:4
    p(:,ii) = desired_pattern(:,angle_test(ii)+51);
end

 w_data_raw = pinv(FF')*sum(p,2);
 

 w_data = w_data_raw./norm(w_data_raw)*sqrt(Nr);


% Pattern Evaluation
backoff_dB = 0;
backoff  = 10^(backoff_dB/10);
sig = (randn(1e3,1) + 1j*randn(1e3,1))/sqrt(2)/sqrt(backoff);

figure
phi = (-180:179)/180*pi;
P = 2;
vsat = 1;
for ss=1:1e3
    sig_NLPA = get_rapp_NL(w_data * sig(ss),1,2);
% response = abs(FF'*(w_data.*exp(1j*(rand(Nr,1)-0.5)*(90/180*pi)))); %with phase error
    response(ss,:) = abs(FF'*sig_NLPA./sqrt(Nr)*sqrt(backoff));
end
response_mean = sqrt(mean(abs(response).^2,1));
plot(phi/pi*180,20*log10(abs(response_mean)));
% polarplot(phi,circshift(response_mean,90));
% polarplot(phi,(abs(FF'*w_data)));

hold on;
% naive_BF_vec = exp(1j*pi*(0:Nr-1).'*sin(angle_test/180*pi))/sqrt(Nr);
% polarplot(phi,20*log10(abs(FF'*naive_BF_vec)));
% legend('codebook','naive')
hold on
grid on
% thetaticks(linspace(0,360-360/36,36))
% thetalim([0,180])
% thetaticklabels({'0','20',})
% ylim([-10,15])