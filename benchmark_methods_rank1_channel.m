% 10/05/2017
% Tracking receiver angle in UE mobility using 
% 1) neighbor angle refinement search
% 2) derivative based narrow band channel tracking
% Channel used here has only one multipath cluster

%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;close all
rng(3); %random seed
print_stat = 1; % print channel parameter summary

%-------------------------------------
% System Parameters
%-------------------------------------
Nt = 32; % Number of Tx antennas (ULA)
Nr = 8; % Number of Tx antennas (ULA)
cluster_num = 1; % Number of multipath clusters
ray_num = 20; % Number of intra-cluster rays
noise_power = 0.1;
sigma_delay_spread = 0;
centroid_AOA = 0/180*pi;
centroid_AOD = 0/180*pi;
sigma_AOA_spread = 5/180*pi;
sigma_AOD_spread = 0/180*pi;

%-------------------------------------
% Generate channel parameters from geolocations of BS, UE, and scatterers
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

%-------------------------------------                         
% Frequency Domain MIMO Channel & Frequency Domain Equivalent Beamspace Channel
%-------------------------------------
H_freq0 = get_H_NB( raygain,...
                    ray_AOA_azim,...
                    ray_AOD_azim,...
                    cluster_num,...
                    ray_num,...
                    Nt, Nr);
                   
%% Narrow Band Compressive Tracking Beacon (D.Ramasamy et el)
% each element is from \pm 1 \pm 1j
M_beacon = 16; % number of beacons
beacon = (randi(2,Nr,M_beacon)*2-3) + 1j*(randi(2,Nr,M_beacon)*2-3);
W = beacon./norm(beacon)*sqrt(Nr*M_beacon);

Probe_precoder = (randi(2,Nt,M_beacon)*2-3)+1j*(randi(2,Nt,M_beacon)*2-3);
for mm = 1:M_beacon
%         Probe_precoder(:,mm) = exp(1j*(0:Nt-1)'*pi*sin(true_rayAOD(1,1))); 
    F(:,mm) = Probe_precoder(:,mm)./norm(Probe_precoder(:,mm))*sqrt(Nt);
end

%% Beam Sweep Based Training
angle_tx_range = linspace(-pi/3,pi/3,Nt*2);
angle_rx_range = linspace(-pi/3,pi/3,Nr*2);
[ tx_max_index,...
  tx_angle,...
  rx_max_index,...
  rx_angle ] = run_angle_sweep_2D( Nt, Nr, 1,...
                                   H_freq0,...
                                   angle_tx_range,...
                                   angle_rx_range);

% figure
% surf(angle_tx_range/pi*180,angle_rx_range/pi*180,rx_RSS)
% xlabel('AOD [deg]')
% ylabel('AOA [deg]')

% surf(1:64,1:32,rx_RSS)

%% Channel change in NB channel
AoA_rotate = 100/180*pi; % 100 deg/s rotation
t_range = (1:1:500)*1e-3;
for tt = 1:length(t_range)
    t_now = t_range(tt);
    dt = t_range(2)-t_range(1);
    
    centroid_AOA(tt+1) = centroid_AOA(tt) + AoA_rotate *dt;
    
    clc
    if tt==1
        raygain(tt,:) = exp(1j*rand(1,ray_num)*2*pi)/sqrt(ray_num);
    else
        raygain(tt,:) = raygain(tt-1,:).*exp(1j*(rand(1,ray_num)*40-20)*(pi/180));
    end
    [ ~, raydelay(tt,:), ray_AOA_azim(tt,:), ray_AOD_azim(tt,:) ] =...
    get_chan_parameter_nogeo( print_stat,...
                              cluster_num,...
                              ray_num,...
                              sigma_delay_spread,...
                              centroid_AOA(tt),...
                              sigma_AOA_spread,...
                              centroid_AOD,...
                              sigma_AOD_spread );
    
end

%% Parameter tracking using beam refinement
% stepsize = 1e-7;
% 
% raydelay_est_BB1 = zeros(length(t_range),ray_num);
% rayAOA_est_BB1 = zeros(length(t_range),ray_num);
% raygain_est_BB1 = zeros(length(t_range),ray_num);

rx_opt_index = zeros(length(t_range)+1,1);
rx_opt_index(1) = rx_max_index;

theta_opt = centroid_AOD;
atx_opt = exp(1j*(0:Nt-1)'*pi*sin(theta_opt))/sqrt(Nt);


angle_est_rad(1) = angle_rx_range(rx_max_index);

for tt = 1:length(t_range)
    
    t_now = t_range(tt);
    H_NB_time_evo = get_H_NB( raygain(tt,:),...
                        ray_AOA_azim(tt,:),...
                        ray_AOD_azim(tt,:),...
                        cluster_num,...
                        ray_num,...
                        Nt, Nr);
%     MIMO_noise = (randn(Nr,Nt,Nfft)+1j*randn(Nr,Nt,Nfft))*sqrt(noise_pow/2);
%     H_freq_noisy = H_freq + MIMO_noise;
    
    %% Neighbor Angle Refinement
    if mod(tt,20)==1
        rx_opt_index(tt+1) =...
            run_angle_refinement_2D( Nr,...
                                     rx_opt_index(tt),...
                                     atx_opt,...
                                     1,...
                                     H_NB_time_evo,...
                                     angle_rx_range);
    else
        rx_opt_index(tt+1) = rx_opt_index(tt);
    end
    
    
    %% NB tracking by updating angle & gain of each cluster
    if mod(tt,20)==1
        noise_vec = randn(M_beacon,1) + 1j*randn(M_beacon,1);
        awgn = noise_vec * sqrt(noise_power/2);
        chan_noisy_ob = diag(W' * H_NB_time_evo * F) + awgn;
        
        angle_est_rad(tt+1) =...
            run_refinement(chan_noisy_ob,...
            angle_est_rad(tt),...
            centroid_AOD,...
            F, W,...
            cluster_num,...
            Nt, Nr);
    else
        angle_est_rad(tt+1) = angle_est_rad(tt);
    end

    
    %% Evaluation of capacity of Neighbor Angle Refinement
    phi_opt = angle_rx_range(rx_opt_index(tt+1));
    arx_opt = exp(1j*(0:Nr-1)'*pi*sin(phi_opt))/sqrt(Nr);
    rx_gain_sector(tt) = arx_opt' * H_NB_time_evo * atx_opt;
    SINR_sector(tt) = abs(rx_gain_sector(tt))^(2) / (noise_power);
    capacity_sector(tt) = log2(1+SINR_sector(tt)); % bps/Hz
    
    %% Evaluation of capacity of NB Tracking
    arx_opt_NBtrack = exp(1j*(0:Nr-1)'*pi*sin(angle_est_rad(tt+1)))/sqrt(Nr);
    rx_gain_NBtrack(tt) = arx_opt_NBtrack' * H_NB_time_evo * atx_opt;
    SINR_NBtrack(tt) = abs(rx_gain_NBtrack(tt))^(2) / (noise_power);
    capacity_NBtrack(tt) = log2(1+SINR_NBtrack(tt)); % bps/Hz
    
    %% Benchmark with perfect CSI

    L = cluster_num;
    [U,S,V] = svd(H_NB_time_evo);
    [diag_S, temp_idx] = sort(diag(S), 'descend');
    S_sorted  = diag(diag_S);
    U_sorted = U(:,temp_idx);
    V_sorted = V(:,temp_idx); 
    U_est = U_sorted(:,1:L);
    V_est = V_sorted(:,1:L);

    % true chan
    gain = U_est' * H_NB_time_evo *V_est;

    if L==2
        % Capcity computing with L = 2 multiplexed streams
        SINR_trueCSI(1, kk) = abs(gain(1,1))^(2) / ((abs(gain(1,2)).^2) + noise_pow);
        SINR_trueCSI(2, kk) = abs(gain(2,2))^(2) / ((abs(gain(2,1)).^2) + noise_pow);
    elseif L==1
        % Capcity computing with L = 2 multiplexed streams
        SINR_trueCSI(tt) = abs(gain)^(2) / (noise_power);
    end

    capacity_trueCSI(tt) = log2(1+SINR_trueCSI(tt)); % bps/Hz


end


%% angle tracking performance
figure
plot(angle_rx_range(rx_opt_index(2:end))/pi*180,'linewidth',2);hold on
plot(angle_est_rad(2:end)/pi*180,'linewidth',2);hold on
plot(centroid_AOA/pi*180,'linewidth',2);hold on
legend('Neighbor Angle Refinement','NB Tracking','True AOA Mean')
grid on

% %%  capacity plotting
% figure;
% plot(t_range, SINR_sector,'--o','linewidth',2);hold on
% plot(t_range, SINR_NBtrack,'--x','linewidth',2);hold on
% plot(t_range, SINR_trueCSI,'linewidth',2);hold on
% legend('Neighbor Angle Refinement','NB Tracking','Perfect CSI')
% grid on
%%  capacity plotting
figure;
plot(t_range, SINR_trueCSI-SINR_sector,'linewidth',2);hold on
plot(t_range, SINR_trueCSI-SINR_NBtrack,'linewidth',2);hold on
legend('Neighbor Angle Refinement','NB Tracking')
title('Gain')
xlabel('Time [s]')
ylabel('Gain Drop [dB]')
grid on
        