% test whether it is possible to get reasonable "initial guess" of channle
% parameter
clear;clc;
ray_num = 20;
MCtimes = 5;
speed_rotate = 100/180*pi;

% Tracking update rate is 250 Hz (one training frame every 4 ms)
t_range = (4:4:500)*1e-3;

% plot_prra control whehher to plot parameter estimation performance of
% each channel realizations
plot_para = zeros(MCtimes,1);
plot_para(1:10) = 1;

capacity_all_OMPstart = zeros(length(t_range), MCtimes);
delay_err_all = zeros(length(t_range), MCtimes*2);
AOA_err_all = zeros(length(t_range), MCtimes*2);

% it is little confusing. Here L indicate number of stream to use.
% Clusyer_num is the L in paper
L = 1;

% SNR_range = [-10,0,10];
SNR_range = [0];
for ss = 1:length(SNR_range)
    noise_pow = 10^(-SNR_range(ss)/10);
    for MCindex=1:MCtimes
        fprintf('MC Realization %d \n',MCindex);

           [capacity_prop, capacity_true, delay_err, AOA_err]...
               = capacity_eval_v3(noise_pow,...
                                  t_range,...
                                  speed_rotate,...
                                  plot_para(MCindex),...
                                  MCindex,...
                                  L);

           capacity_all_OMPstart(:,MCindex) = capacity_prop;
           capacity_all_true(:,MCindex) = capacity_true;

           delay_err_all(:,(MCindex-1)*2+1:MCindex*2) = delay_err;
           AOA_err_all(:,(MCindex-1)*2+1:MCindex*2) = AOA_err;
    end
    
    %  capacity plotting
    figure
    plot(t_range, mean(capacity_all_OMPstart,2));hold on
    plot(t_range, mean(capacity_all_true,2));hold on
    legend('SVD from Est.','SVD from True Chan.')
    xlabel('Time (second)')
    ylabel('Spactral Efficiency (bps/hz)')
    titlename = ['Capcity, SNR = ',num2str(SNR_range(ss)),' dB'];
    title(titlename)

    % parameter plotting
    delay_err_rms = sqrt(mean(abs(delay_err_all).^2,2))*1e9;
    AOA_err_rms = sqrt(mean(abs(AOA_err_all).^2,2))/pi*180;
    
    figure
    subplot(211)
    plot(t_range, delay_err_rms);hold on
    xlabel('Time (second)')
    ylabel('Delay (nanosecond)')
    titlename = ['RMS Delay Est. Error, SNR = ',num2str(SNR_range(ss)),' dB'];
    title(titlename)

    subplot(212)
    plot(t_range, AOA_err_rms);
    xlabel('Time (second)')
    ylabel('Angle (degree)')
    titlename = ['RMS AOA Est. Error, SNR = ',num2str(SNR_range(ss)),' dB'];
    title(titlename)
end

%%
figure
plot(t_range, mean(capacity_all_OMPstart(:,1:5),2));hold on
plot(t_range, mean(capacity_all_true,2));hold on