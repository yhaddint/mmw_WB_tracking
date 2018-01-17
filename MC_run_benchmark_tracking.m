% test whether it is possible to get reasonable "initial guess" of channle
% parameter
clear;clc;
ray_num = 20;
MCtimes = 5;
speed_rotate = 50/180*pi;

% Tracking update rate is 250 Hz (one training frame every 4 ms)
t_range = (10:10:500)*1e-3;

% plot_prra control whehher to plot parameter estimation performance of
% each channel realizations
plot_para = zeros(MCtimes,1);
plot_para(1:10) = 1;

capacity_all_refine = zeros(length(t_range), MCtimes);
capacity_all_NBtrack = zeros(length(t_range), MCtimes);
capacity_all_trueCSI = zeros(length(t_range), MCtimes);

% it is little confusing. Here L indicate number of stream to use.
% Clusyer_num is the L in paper
L = 1;

% SNR_range = [-10,0,10];
SNR_range = [0];
Nr = 16;
M_beacon = 30; % number of beacons
beacon = ((randi(2,Nr,M_beacon)*2-3) + 1j*(randi(2,Nr,M_beacon)*2-3))...
        /sqrt(Nr*2);
    
for ss = 1:length(SNR_range)
    noise_pow = 10^(-SNR_range(ss)/10);
    for MCindex=1:MCtimes
        fprintf('MC Realization %d \n',MCindex);

           [capacity_refine, capacity_NBtrack, capacity_trueCSI]...
               = capacity_benchmark_v1(noise_pow,...
                                       t_range,...
                                       speed_rotate,...
                                       plot_para(MCindex),...
                                       MCindex,...
                                       L,...
                                       beacon);

           capacity_all_refine(:,MCindex) = capacity_refine;
           capacity_all_NBtrack(:,MCindex) = capacity_NBtrack;
           capacity_all_trueCSI(:,MCindex) = capacity_trueCSI;

%            delay_err_all(:,(MCindex-1)*2+1:MCindex*2) = delay_err;
%            AOA_err_all(:,(MCindex-1)*2+1:MCindex*2) = AOA_err;
    end
    
    %  capacity plotting
    figure
    plot(t_range, mean(capacity_all_refine,2));hold on
    plot(t_range, mean(capacity_all_NBtrack,2));hold on
    plot(t_range, mean(capacity_all_trueCSI,2));hold on
    legend('Refine','NB track','True CSI')
    xlabel('Time (second)')
    ylabel('Spactral Efficiency (bps/hz)')
    titlename = ['Capcity, SNR = ',num2str(SNR_range(ss)),' dB'];
    title(titlename)

    % parameter plotting
%     delay_err_rms = sqrt(mean(abs(delay_err_all).^2,2))*1e9;
%     AOA_err_rms = sqrt(mean(abs(AOA_err_all).^2,2))/pi*180;
    
%     figure
%     subplot(211)
%     plot(t_range, delay_err_rms);hold on
%     xlabel('Time (second)')
%     ylabel('Delay (nanosecond)')
%     titlename = ['RMS Delay Est. Error, SNR = ',num2str(SNR_range(ss)),' dB'];
%     title(titlename)
% 
%     subplot(212)
%     plot(t_range, AOA_err_rms);
%     xlabel('Time (second)')
%     ylabel('Angle (degree)')
%     titlename = ['RMS AOA Est. Error, SNR = ',num2str(SNR_range(ss)),' dB'];
%     title(titlename)
end

%%
figure
plot(t_range, mean(capacity_all_refine,2));hold on
plot(t_range, mean(capacity_all_NBtrack(:,[1:5]),2));hold on
plot(t_range, mean(capacity_all_trueCSI,2));hold on