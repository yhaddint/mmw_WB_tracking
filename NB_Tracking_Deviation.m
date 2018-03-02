%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;
rng(3)
plot_ellipse = 0;
print_stat = 1;

%-------------------------------------
% System Parameters
%-------------------------------------
Nt = 32; % Number of Tx antennas (ULA)
Nr = 8; % Number of Tx antennas (ULA)
cluster_num = 1; % Number of multipath clusters
ray_num = 10; % Number of intra-cluster rays
sigma_delay_spread = 0;
centroid_AOA = 0/180*pi;
sigma_AOA_spread = 0/180*pi;
centroid_AOD = 0;
sigma_AOD_spread = 0;
M = 16;
MCtimes = 200;

%% ML estimation of angle

noise_pow_range = -15:2.5:30;
errorSD_noisy = zeros(length(noise_pow_range),1); %RMSE

for MCindex = 1:MCtimes %1000 Monte Carlo simulations
    
    % Generate channel parameters from desired statistics
    clc
    fprintf('Monte Carlo ite %4.0f\n',MCindex);
    [ raygain, raydelay, ray_AOA_azim, ray_AOD_azim ] =...
        get_chan_parameter_nogeo( print_stat,...
                                  cluster_num,...
                                  ray_num,...
                                  sigma_delay_spread,...
                                  centroid_AOA,...
                                  sigma_AOA_spread,...
                                  centroid_AOD,...
                                  sigma_AOD_spread );
    
    % Frequency Domain MIMO Channel & Frequency Domain Equivalent Beamspace Channel
    H_NB = get_H_NB( raygain,...
                     ray_AOA_azim,...
                     ray_AOD_azim,...
                     cluster_num,...
                     ray_num,...
                     Nt, Nr);
    H0 = H_NB;
%     norm_factor = norm(H0,'fro');
    H_normal = H0;%/norm_factor*sqrt(Nt*Nr);

    
    sig = 0; %received signal
    true_rayAOA = zeros(1, cluster_num);
    true_rayAOD = zeros(1, cluster_num);
        
    Probe_combiner = (randi(2,Nr,M)*2-3)+1j*(randi(2,Nr,M)*2-3);
    for mm = 1:M
        W(:,mm) = Probe_combiner(:,mm)./norm(Probe_combiner(:,mm))*sqrt(Nr);
    end
    
%     Probe_precoder = (randi(2,Nt,M)*2-3)+1j*(randi(2,Nt,M)*2-3);
    for mm = 1:M
        Probe_precoder(:,mm) = exp(1j*(0:Nt-1)'*pi*sin(true_rayAOD(1,1))); 
        F(:,mm) = Probe_precoder(:,mm)./norm(Probe_precoder(:,mm))*sqrt(Nt);
    end
        
        
    for cluster_index = 1:cluster_num
        true_rayAOA(1,cluster_index) = centroid_AOA;
        true_rayAOD(1,cluster_index) = centroid_AOD;
%         Probe_predocer = exp(1j*(0:Nt-1)'*pi*sin(true_rayAOD(1,cluster_index))); %we assume transmitter is perfectly aligned

        %if we have multiple clusters, must save information about true
        %rayAOA, true rayAOD and steering vectors (e.g. make new vectors
        %and save them into these vectors)
%         F = repmat(Probe_predocer, 1, M); %we assume that receiver makes 4 measurements (for example)
        sig = sig + diag(W'*H_normal*F); %we take just one subcarrier
    end
        
    % For loop for different SNR
    noise_vec = randn(M,1) + 1j*randn(M,1);
    angle_est_rad = zeros(length(noise_pow_range), 1);
    for noise_index = 1:length(noise_pow_range)
        noise_power = 10^(-noise_pow_range(noise_index)/10);
        awgn = noise_vec * sqrt(noise_power/2);
        chan_noisy_ob = sig + awgn;
        angle_est_rad(noise_index) = ml_angle(chan_noisy_ob, true_rayAOA, true_rayAOD, F, W, cluster_num, Nt, Nr);
        
%         % Evaluate SNR gain
%         [U,Sigma,V] = svd(H_normal);
%         data_combiner = exp(1j*pi*(0:Nr-1).'*sin(angle_est_rad(noise_index)));
%         data_combiner_theo = exp(1j*pi*(0:Nr-1).'*sin(centroid_AOA(1,1)+randn*(5/180*pi)));
%         data_combiner_opt = U(:,1)./abs(U(:,1));
%         data_precoder = exp(1j*pi*(0:Nt-1).'*sin(true_rayAOD(1,1)));
%         gain_ML(noise_index,MCindex) = abs(data_combiner'*H_normal*data_precoder)^2;
%         gain_theo(noise_index,MCindex) = abs(data_combiner_theo'*H_normal*data_precoder)^2;
%         gain_opt(noise_index,MCindex) = abs(data_combiner_opt'*H_normal*data_precoder)^2;
    end
    errorSD_noisy(:,MCindex) = (angle_est_rad - true_rayAOA(1,1)).^2;
    %angle_error_var = var(angle_est);
end
%% Figure plotting
MSE_noisy = sqrt(mean(errorSD_noisy,2))/pi*180;
figure
semilogy(noise_pow_range, MSE_noisy, 'Linewidth', 2);
xlabel('SNR (dB)')
ylabel('RMSE of AOA (deg)')
title('RMSE vs. SNR[dB]')
grid on

%%
% gain_opt_mean = mean(gain_opt,2);
% gain_mean = mean(gain_ML,2);
% gain_theo_mean = mean(gain_theo,2);
% figure
% plot(noise_pow_range,10*log10(gain_opt_mean));hold on
% plot(noise_pow_range,10*log10(gain_mean));hold on
% plot(noise_pow_range,10*log10(gain_theo_mean));hold on
% 
% grid on
% legend('true centroid','ML','from CRLB')
