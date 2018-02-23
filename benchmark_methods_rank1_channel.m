% 10/05/2017
% Tracking receiver angle in UE mobility using 
% 1) neighbor angle refinement search
% 2) derivative based narrow band channel tracking
% Channel used here has only one multipath cluster

%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;close all
rng(2); %random seed
plot_ellipse = 0; % plot geolocations of Tx/Rx/Scatterers
print_stat = 0; % print channel parameter summary

%-------------------------------------
% System Parameters
%-------------------------------------
Nt = 32; % Number of Tx antennas (ULA)
Nr = 16; % Number of Tx antennas (ULA)
fc = 28e9; % Carrier frequency
cluster_num = 1; % Number of multipath clusters
ray_num = 20; % Number of intra-cluster rays
noise_pow = 10;
Nfft = 512; % Number of subcarriers

%-------------------------------------
% Generate geolocations of BS, UE, and scatterers
%-------------------------------------
[ loc0_ue,...
  loc0_bs,...
  loc_cluster_total ] = get_init_locations( plot_ellipse,...
                                            cluster_num,...
                                            ray_num);

%-------------------------------------
% Generate channel parameters from geolocations of BS, UE, and scatterers
%-------------------------------------
[raydelay,...
 rayAOA,...
 rayAOD ] = get_multipath_v2(loc0_bs,...
                             loc0_ue,...
                             loc_cluster_total,...
                             cluster_num,...
                             ray_num,...
                             print_stat);
raygain = exp(1j*rand(cluster_num, ray_num)*2*pi); % Complex gain is randomly generated

%-------------------------------------                         
% Frequency Domain MIMO Channel & Frequency Domain Equivalent Beamspace Channel
%-------------------------------------
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

%% Narrow Band Compressive Tracking Beacon (D.Ramasamy et el)
% each element is from \pm 1 \pm 1j
M_beacon = 30; % number of beacons
beacon = ((randi(2,Nr,M_beacon)*2-3) + 1j*(randi(2,Nr,M_beacon)*2-3))...
        /sqrt(Nr*2);

%% Beam Sweep Based Training
angle_tx_range = linspace(-pi/2,pi/2,Nt*2);
angle_rx_range = linspace(-pi/2,pi/2,Nr*2);
[ tx_max_index,...
  tx_angle,...
  rx_max_index,...
  rx_angle ] = run_angle_sweep_2D( Nt, Nr, Nfft,...
                                   H_freq0,...
                                   angle_tx_range,...
                                   angle_rx_range);

% figure
% surf(angle_tx_range/pi*180,angle_rx_range/pi*180,rx_RSS)
% xlabel('AOD [deg]')
% ylabel('AOA [deg]')

% surf(1:64,1:32,rx_RSS)

%% Channel change in BB1 and BB2
rho = 0.999;
speed_rotate = -100/180*pi;
t_range = (4:4:200)*1e-3;
for tt = 1:length(t_range)
    t_now = t_range(tt);
    dt = t_range(2)-t_range(1);
    
    % setup speed
    if t_now<500e-3
        speed_v = [10, 0];
    else
        speed_v = [10, 0];
    end
    
    if tt==1
        loc_ue = loc0_ue + speed_v * dt;
    else
        loc_ue = loc_ue + speed_v * dt;
    end
     
%     clc
%     fprintf('Time Evolution %2.4f s\n',t_range(tt));
    [raydelay,...
     rayAOA,...
     rayAOD ] = get_multipath_v2(loc0_bs,...
                                 loc_ue,...
                                 loc_cluster_total,...
                                 cluster_num,...
                                 ray_num,...
                                 print_stat);
                             
    rayAOA = rayAOA + speed_rotate * (tt-1) * dt;
    raygain = raygain.*exp(1j*rand(cluster_num,ray_num)*2*pi*sqrt(1-rho^2));

    raygain_true_BB1(tt,:)  = raygain(1,:);
    raydelay_true_BB1(tt,:)  = raydelay(1,:);
    rayAOA_true_BB1(tt,:)  = rayAOA(1,:);
end

%% Parameter tracking using beam refinement
% stepsize = 1e-7;
% 
% raydelay_est_BB1 = zeros(length(t_range),ray_num);
% rayAOA_est_BB1 = zeros(length(t_range),ray_num);
% raygain_est_BB1 = zeros(length(t_range),ray_num);

rx_opt_index = zeros(length(t_range)+1,1);
rx_opt_index(1) = rx_max_index;

theta_opt = angle_tx_range(tx_max_index);
atx_opt = exp(1j*(0:Nt-1)'*pi*sin(theta_opt))/sqrt(Nt);


for tt = 1:length(t_range)
    
    raygain(1,:) = raygain_true_BB1(tt,:);
    raydelay(1,:) = raydelay_true_BB1(tt,:);
    rayAOA(1,:) = rayAOA_true_BB1(tt,:);
    
    H_freq = get_H_freq2(raygain,...
                         raydelay,...
                         rayAOA,...
                         rayAOD,...
                         cluster_num,...
                         ray_num,...
                         Nt, Nr);
    H_freq = H_freq / norm_factor;
%     MIMO_noise = (randn(Nr,Nt,Nfft)+1j*randn(Nr,Nt,Nfft))*sqrt(noise_pow/2);
%     H_freq_noisy = H_freq + MIMO_noise;
    
    %% Neighbor Angle Refinement
    rx_opt_index(tt+1) = ...
        run_angle_refinement_2D( Nr,...,
                                 rx_opt_index(tt),...
                                 atx_opt,...
                                 Nfft,...
                                 H_freq,...
                                 angle_rx_range);
    
    
    %% NB tracking by updating angle & gain of each cluster
    %---------------------------------------------
    % Alpha estimation using previous Phi
    %---------------------------------------------
    if tt==1
        phi_hat = mean(rayAOA(1,:),2); % Assume mean angle is known (practically via CS)
    else
        phi_hat = phi_new;
    end
    kk = 1;
    rx_measure = zeros(M_beacon,1);
    Hk = squeeze(H_freq(:,:,kk));
    for m_beacon = 1:M_beacon
        rx_measure(m_beacon) = beacon(:,m_beacon)' * Hk * atx_opt;
    end
    arx_hat = exp(1j*(0:Nr-1)'*pi*sin(phi_hat))/sqrt(Nr);
    H_cal = beacon' * arx_hat;
    alpha_est = pinv(H_cal) * rx_measure;
    
    % quick watch received beacon
%     figure;
%     plot(abs(H_cal * alpha_est));hold on
%     plot(abs(rx_measure))
    
    %---------------------------------------------
    % Phi estimation using previous Alpha
    %---------------------------------------------
    dPhi = zeros(Nr,1);
    for nn=1:Nr
        dPhi(nn) = exp(1j*pi*(nn-1)*sin(phi_hat))/sqrt(Nr)*(1j*pi*(nn-1)*cos(phi_hat));
    end
    D_cal = alpha_est * (beacon' * dPhi);
    yr = rx_measure - H_cal * alpha_est;
    Q_cal = [real(D_cal);imag(D_cal)];
    tr = [real(yr);imag(yr)];
    
    dphi = pinv(Q_cal) * tr;
    phi_new = phi_hat + dphi;
    phi_array(tt) = phi_new;
    
    %% Evaluation of capacity of Neighbor Angle Refinement
    phi_opt = angle_rx_range(rx_opt_index(tt+1));
    arx_opt = exp(1j*(0:Nr-1)'*pi*sin(phi_opt))/sqrt(Nr);
    for kk=1:Nfft
        H_k = squeeze(H_freq(:,:,kk));
        rx_gain(kk) = arx_opt' * H_k * atx_opt;
        SINR_refine(kk) = abs(rx_gain(kk))^(2) / (noise_pow);
    end
    capacity_refine(tt) = mean(log2(1+SINR_refine)); % bps/Hz
    
    %% Evaluation of capacity of NB Tracking
    arx_opt_NBtrack = exp(1j*(0:Nr-1)'*pi*sin(phi_new))/sqrt(Nr);
    for kk=1:Nfft
        H_k = squeeze(H_freq(:,:,kk));
        rx_gain(kk) = arx_opt_NBtrack' * H_k * atx_opt;
        SINR_NBtrack(kk) = abs(rx_gain(kk))^(2) / (noise_pow);
    end
    capacity_NBtrack(tt) = mean(log2(1+SINR_NBtrack)); % bps/Hz

    %% Benchmark with perfect CSI
    for kk=1:Nfft
        L = cluster_num;
        [U,S,V] = svd(H_freq(:,:,kk));
        [diag_S, temp_idx] = sort(diag(S), 'descend');
        S_sorted  = diag(diag_S);
        U_sorted = U(:,temp_idx);
        V_sorted = V(:,temp_idx); 
        U_est{kk} = U_sorted(:,1:L);
        V_est{kk} = V_sorted(:,1:L);
        
        % true chan
        H_k = squeeze(H_freq(:,:,kk));
        gain = U_est{kk}' * H_k *V_est{kk};
        
        if L==2
            % Capcity computing with L = 2 multiplexed streams
            SINR_trueCSI(1, kk) = abs(gain(1,1))^(2) / ((abs(gain(1,2)).^2) + noise_pow);
            SINR_trueCSI(2, kk) = abs(gain(2,2))^(2) / ((abs(gain(2,1)).^2) + noise_pow);
        elseif L==1
            % Capcity computing with L = 2 multiplexed streams
            SINR_trueCSI(1, kk) = abs(gain(1,1))^(2) / (noise_pow);
        end
    end
    if L==2
        capacity_trueCSI(tt) = mean(sum(log2(1+SINR_trueCSI))); % bps/Hz
    elseif L==1
        capacity_trueCSI(tt) = mean(log2(1+SINR_trueCSI)); % bps/Hz
    end

end


%% angle tracking performance
figure
plot(angle_rx_range(rx_opt_index)/pi*180,'linewidth',2);hold on
plot(phi_array/pi*180,'linewidth',2);hold on
plot(mean(rayAOA_true_BB1,2)/pi*180,'linewidth',2);hold on
legend('Neighbor Angle Refinement','NB Tracking','True AOA Mean')
grid on

%%  capacity plotting
figure;
plot(t_range, capacity_refine,'linewidth',2);hold on
plot(t_range, capacity_NBtrack,'linewidth',2);hold on
plot(t_range, capacity_trueCSI,'linewidth',2);hold on
legend('Neighbor Angle Refinement','NB Tracking','Perfect CSI')
grid on
        