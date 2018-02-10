%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc
rng(2)
plot_ellipse = 0;
print_stat = 1;

%-------------------------------------
% System Parameters
%-------------------------------------
Nt = 32; % Number of Tx antennas (ULA)
Nr = 32; % Number of Tx antennas (ULA)
fc = 28e9; % Carrier frequency
cluster_num = 1; % Number of multipath clusters
ray_num = 10; % Number of intra-cluster rays
sigma_delay_spread = 0;
centroid_AOA = 30/180*pi;
sigma_AOA_spread = 5/180*pi;
centroid_AOD = 0;
sigma_AOD_spread = 0;
Nfft = 512; % Number of subcarriers
M = 4;
MCtimes = 100;
%-------------------------------------
% Generate channel parameters from desired statistics
%-------------------------------------
[ raygain, raydelay, ray_AOA_azim, ray_AOD_azim ] =...
    get_chan_parameter_nogeo( print_stat,...
                              cluster_num,...
                              ray_num,...
                              sigma_delay_spread,...
                              centroid_AOA,...
                              sigma_AOA_spread,...
                              centroid_AOD,...
                              sigma_AOD_spread );
                          

%% Frequency Domain MIMO Channel & Frequency Domain Equivalent Beamspace Channel
Nfft = 512;
H_freq0 = get_H_freq2( raygain,...
                       raydelay,...
                       ray_AOA_azim,...
                       ray_AOD_azim,...
                       cluster_num,...
                       ray_num,...
                       Nt, Nr);
norm_factor = sqrt(mean(mean(mean(abs(H_freq0).^2))));
% norm_factor = 1;
H_freq0 = H_freq0 / norm_factor ;

%% ML estimation of angle
%AOA = rand*pi-pi/2;
noises = -20:5:20;
errorSD_noisy = zeros(length(noises),1); %RMSE
for noise_index = 1:length(noises);
    angle_est = zeros(1000, 1);
    for MCindex = 1:MCtimes %1000 Monte Carlo simulations
        y = 0; %received signal
        true_rayAOA = zeros(1, cluster_num);
        true_rayAOD = zeros(1, cluster_num);
        for cluster_index = 1:cluster_num
            true_rayAOA(1,cluster_index) = centroid_AOA;
            true_rayAOD(1,cluster_index) = centroid_AOD;
            Probe_predocer = exp(1j*(0:Nt-1)'*pi*sin(true_rayAOD(1,cluster_index)))/sqrt(Nt); %we assume transmitter is perfectly aligned

            %if we have multiple clusters, must save information about true
            %rayAOA, true rayAOD and steering vectors (e.g. make new vectors
            %and save them into these vectors)
            F = repmat(Probe_predocer, 1, M); %we assume that receiver makes 4 measurements (for example)
            Probe_combiner = (randi(2,Nr,4)*2-3)+1j*(randi(2,Nr,4)*2-3);
            W = Probe_combiner./norm(Probe_combiner,'fro');
            y = y + W'*H_freq0(:,:,1)*F*ones(4,1); %we take just one subcarrier
            signal_power = norm(y)^2;
            noise_power = signal_power / 10^(noises(noise_index)/10);
            noise_vec = (randn(M,1)*sqrt(noise_power/8) + 1i*randn(M,1)*sqrt(noise_power/8));
            y = y + noise_vec;
        end
        angle_est(MCindex,1) = ml_angle(y, true_rayAOA, true_rayAOD, F, W, cluster_num, Nt, Nr);
    end
    errorSD_noisy(noise_index,1) = sqrt(mean((angle_est - ones(1000,1)*(true_rayAOA(1,1)/pi*180)).^2));
    %angle_error_var = var(angle_est);
end

close all
figure
semilogy(noises, errorSD_noisy, 'o-r', 'Linewidth', 2)
xlabel('SNR [dB]')
ylabel('RMSE of AOA estimation')
title('RMSE vs. SNR[dB]')


