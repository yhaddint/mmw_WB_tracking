% performance with perfect initial parameter...
clear;clc;close all
% rng(2)
%% other parameter
print_stat = 0;

%% Goemetry Parameters
Nt = 16;
Nr = 16;
D_bs2ue = 40;
loc0_ue = [-D_bs2ue/2,0];
loc0_bs = [D_bs2ue/2,0];
speed_c = 3e8;
%% Large Scale Parameter (Multiple Clusters)
cluster_num = 1;
ray_num = 2;
sigma_AOA_spread = 2:2:18; %different deviations for AOA
SNR = 20; %in dB
mc = 1000;

errorSD_noisy = zeros(length(sigma_AOA_spread),1); %RMSE
for spread_index = 1:length(sigma_AOA_spread)
    %% Channel parameters
    [raygain, raydelay, rayAOA, rayAOD, true_rayAOA, true_rayAOD] = get_chan_parameter_nogeo( print_stat,...
                                                                                              cluster_num,...
                                                                                              ray_num,...
                                                                                              0,...
                                                                                              0,...
                                                                                              sigma_AOA_spread(spread_index),...
                                                                                              0,...
                                                                                              0 );

    %% Frequency Domain MIMO Channel & Frequency Domain Equivalent Beamspace Channel
    Nfft = 512;
    H_freq0 = get_H_freq2( raygain,...
                           raydelay,...
                           rayAOA,...
                           rayAOD,...
                           cluster_num,...
                           ray_num,...
                           Nt, Nr);
    norm_factor = sqrt(mean(mean(mean(abs(H_freq0).^2))));
    % norm_factor = 1;
    H_freq0 = H_freq0 / norm_factor ;

    %% ML estimation of angle
    %AOA = rand*pi-pi/2;
    for noise_index = 1:length(SNR);
        angle_est = zeros(mc, 1);
        for i = 1:mc %1000 Monte Carlo simulations
            y = 0; %received signal
            for cluster_index = 1:cluster_num
                atx = exp(1j*(0:Nt-1)'*pi*sin(true_rayAOD(1,cluster_index)))/sqrt(Nt); %we assume transmitter is perfectly aligned
                %arx = exp(1j*(0:Nr-1)'*pi*sin(true_rayAOA(1,cluster_index)))/sqrt(Nr); %we assume transmitter is perfectly aligned
                F = repmat(atx,1,8);%[atx atx atx atx]; %we assume that receiver makes 4 measurements (for example)
                W = (randi(2,Nr,8)*2-3)+1j*(randi(2,Nr,8)*2-3);
                y = y + W'*H_freq0(:,:,1)*F*ones(8,1); %we take just one subcarrier
                signal_power = norm(y)^2;
                noise_power = signal_power / 10^(SNR(noise_index)/10);
                noise_vec = (randn(8,1)*sqrt(noise_power/16) + 1i*randn(8,1)*sqrt(noise_power/16));
                y = y + noise_vec;
            end
            angle_est(i,1) = ml_angle(y, true_rayAOA, true_rayAOD, F, W, cluster_num, Nt, Nr);
        end
        errorSD_noisy(spread_index,1) = sqrt(mean((angle_est - ones(mc,1)*(true_rayAOA(1,1)/pi*180)).^2));
        %angle_error_var = var(angle_est);
    end
end
figure
semilogy(sigma_AOA_spread, errorSD_noisy, '-', 'Linewidth', 1.5)
xlabel('AOA spread')
ylabel('RMSE of AOA estimation')
title('RMSE vs. AOA spread')