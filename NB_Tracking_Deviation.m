%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;
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
ray_num = 4; % Number of intra-cluster rays
sigma_delay_spread = 0;
centroid_AOA = 30/180*pi;
sigma_AOA_spread = 2.5/180*pi;
centroid_AOD = 0;
sigma_AOD_spread = 0;
Nfft = 512; % Number of subcarriers
M = 4;
MCtimes = 100;


%% ML estimation of angle

noise_pow_range = -30:2.5:0;
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
    H_freq0 = get_H_freq2( raygain,...
                           raydelay,...
                           ray_AOA_azim,...
                           ray_AOD_azim,...
                           cluster_num,...
                           ray_num,...
                           Nt, Nr);
    H0 = squeeze(H_freq0(:,:,1));
    norm_factor = norm(H0,'fro');
    H_normal = H0/norm_factor*sqrt(Nt*Nr);


    sig = 0; %received signal
    true_rayAOA = zeros(1, cluster_num);
    true_rayAOD = zeros(1, cluster_num);
        
    Probe_combiner = (randi(2,Nr,M)*2-3)+1j*(randi(2,Nr,M)*2-3);
    for mm = 1:M
        W(:,mm) = Probe_combiner(:,mm)./norm(Probe_combiner(:,mm))*sqrt(Nr);
    end
        
        
    for cluster_index = 1:cluster_num
        true_rayAOA(1,cluster_index) = centroid_AOA;
        true_rayAOD(1,cluster_index) = centroid_AOD;
        Probe_predocer = exp(1j*(0:Nt-1)'*pi*sin(true_rayAOD(1,cluster_index)))/sqrt(Nt); %we assume transmitter is perfectly aligned

        %if we have multiple clusters, must save information about true
        %rayAOA, true rayAOD and steering vectors (e.g. make new vectors
        %and save them into these vectors)
        F = repmat(Probe_predocer, 1, M); %we assume that receiver makes 4 measurements (for example)
        sig = sig + W'*H_normal*F*ones(4,1); %we take just one subcarrier
    end
        
    % For loop for different SNR
    noise_vec = randn(Nr,M) + 1j*randn(Nr,M);
    angle_est_rad = zeros(length(noise_pow_range), 1);
    for noise_index = 1:length(noise_pow_range)
        noise_power = 10^(-noise_pow_range(noise_index)/10);
        awgn = diag(W'*noise_vec*sqrt(noise_power/2));
        chan_noisy_ob = sig + awgn;
        angle_est_rad(noise_index) = ml_angle(chan_noisy_ob, true_rayAOA, true_rayAOD, F, W, cluster_num, Nt, Nr);
    end
    errorSD_noisy(:,MCindex) = (angle_est_rad - true_rayAOA(1,1)).^2;
    %angle_error_var = var(angle_est);
end
%% Figure plotting
MSE_noisy = sqrt(mean(errorSD_noisy,2))/pi*180;
figure
semilogy(noise_pow_range+15.05, MSE_noisy, 'Linewidth', 2);
xlabel('SNR (dB)')
ylabel('RMSE of AOA (deg)')
title('RMSE vs. SNR[dB]')
grid on

