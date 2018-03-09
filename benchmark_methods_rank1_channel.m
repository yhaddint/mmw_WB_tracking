% 10/05/2017
% Tracking receiver angle in UE mobility using 
% 1) neighbor angle refinement search
% 2) derivative based narrow band channel tracking
% Channel used here has only one multipath cluster

%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;
rng(4); %random seed
print_stat = 1; % print channel parameter summary
plot_ellipse = 0;
analyze_angular_profile = 1;
%-------------------------------------
% System Parameters
%-------------------------------------
Nt = 32; % Number of Tx antennas (ULA)
Nr = 16; % Number of Tx antennas (ULA)
fc = 28e9;
cluster_num = 1; % Number of multipath clusters
ray_num = 20; % Number of intra-cluster rays
noise_power = 1;
scatter_radius = 0;
% sigma_delay_spread = 0;
centroid_AOA = 0/180*pi;
centroid_AOD = 0/180*pi;
% sigma_AOA_spread = 5/180*pi;
% sigma_AOD_spread = 0/180*pi;

%-------------------------------------
% Generate geolocations of BS, UE, and scatterers
%-------------------------------------
[ loc0_ue,...
  loc0_bs,...
  loc_cluster_total ] = get_init_locations( plot_ellipse,...
                                            cluster_num,...
                                            ray_num,0);
                                           
%-------------------------------------
% Generate channel parameters from geolocations of BS, UE, and scatterers
%-------------------------------------
[raydelay,...
 ray_AOA_azim,...
 ray_AOD_azim ] = get_multipath_v2(loc0_bs,...
                             loc0_ue,...
                             loc_cluster_total,...
                             cluster_num,...
                             ray_num,...
                             print_stat);
raygain = exp(1j*raydelay*fc*2*pi)/sqrt(ray_num); % Complex gain is randomly generated

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
angle_rx_range = linspace(-pi/3,pi/3,Nr);
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
%% codebook
beam_width = 3;
angle_range = 80;
pointing_range = -50:1:50;
desired_pattern = zeros(angle_range*2+1,length(pointing_range));
for ii=1:length(pointing_range)
    left_boundary = pointing_range(ii) - (beam_width-1)/2;
    right_boundary = pointing_range(ii) + (beam_width-1)/2;
    index_desired = (left_boundary:right_boundary) + angle_range + 1;
    desired_pattern(index_desired,ii) = sqrt(Nr)*ones(beam_width,1);
end
for kk = 1:(angle_range*2+1)
    FF(:,kk) = exp(1j*pi*(0:Nr-1).'*sin((kk - angle_range -1 )/180*pi));
end

for ii=1:101
     w_data_raw = pinv(FF')*desired_pattern(:,ii);
     w_data(:,ii) = w_data_raw./norm(w_data_raw);
end
%%
% angle_test = 0;
% figure
% plot(-80:80,20*log10(abs(FF'*w_data(:,angle_test+51))));
% hold on;
% naive_BF_vec = exp(1j*pi*(0:Nr-1).'*sin(angle_test/180*pi))/sqrt(Nr);
% plot(-80:80,20*log10(abs(FF'*naive_BF_vec)));
% legend('codebook','naive')
% hold on
% grid on
% ylim([0,15])
%% Channel change in NB channel
cluster_AOA(1) = 0;
t_range = (1:1:1000)*1e-3;
for tt = 1:length(t_range)
    t_now = t_range(tt);
    dt = t_range(2)-t_range(1);
    
 % setup speed
    if t_now<500e-3
        speed_v = [5, 0];
        speed_rotate = -50/180*pi;
    else
        speed_v = [0, 5];
        speed_rotate = 25/180*pi;
    end
    
    if tt==1
        loc_ue = loc0_ue + speed_v * dt;
    else
        loc_ue = loc_ue + speed_v * dt;
    end
     
%     clc
%     fprintf('Time Evolution %2.4f s\n',t_range(tt));
    [raydelay,...
     ray_AOA_azim(tt,:),...
     ray_AOD_azim ] = get_multipath_v2(loc0_bs,...
                                 loc_ue,...
                                 loc_cluster_total,...
                                 cluster_num,...
                                 ray_num,...
                                 print_stat);
    
    cluster_AOA(tt+1) = cluster_AOA(tt) + speed_rotate * dt;
    if scatter_radius==0
        ray_AOA_azim(tt,:) = ones(1,ray_num)*mean(ray_AOA_azim(tt,:)) + cluster_AOA(tt+1);
    else
        ray_AOA_azim(tt,:) = ray_AOA_azim(tt,:) + cluster_AOA(tt+1);
    end
%     raygain(tt,:) = exp(1j*(raydelay*fc*2*pi+rand(1,ray_num)*2*pi))/sqrt(ray_num); % Complex gain is randomly generated
    raygain(tt,:) = 1/ray_num; % LOS gain has no fading 

end
%% Power Angular Profile


%% Parameter tracking using beam refinement
% stepsize = 1e-7;
% 
% raydelay_est_BB1 = zeros(length(t_range),ray_num);
% rayAOA_est_BB1 = zeros(length(t_range),ray_num);
% raygain_est_BB1 = zeros(length(t_range),ray_num);

rx_opt_index = zeros(length(t_range)+1,1);
rx_opt_index(1) = rx_max_index;

theta_opt = mean(ray_AOD_azim);
atx_opt = exp(1j*(0:Nt-1)'*pi*sin(theta_opt));


angle_est_rad(1) = angle_rx_range(rx_max_index);

for tt = 1:length(t_range)
    
    t_now = t_range(tt);
    H_NB_time_evo = get_H_NB( raygain(tt,:),...
                        ray_AOA_azim(tt,:),...
                        ray_AOD_azim,...
                        cluster_num,...
                        ray_num,...
                        Nt, Nr);
%     MIMO_noise = (randn(Nr,Nt,Nfft)+1j*randn(Nr,Nt,Nfft))*sqrt(noise_pow/2);
%     H_freq_noisy = H_freq + MIMO_noise;
    
    %%
    if analyze_angular_profile
        power_AOA(:,tt) = analyze_AOA_spead( H_NB_time_evo, Nt, Nr);
    end
    
    %% Neighbor Angle Refinement
    if mod(tt,20)==1
        rx_opt_index(tt+1) =...
            run_angle_refinement_2D( Nr,...
                                     rx_opt_index(tt),...
                                     atx_opt,...
                                     1,...
                                     H_NB_time_evo,...
                                     angle_rx_range,...
                                     0);
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
%     arx_opt = exp(1j*(0:Nr-1)'*pi*sin(phi_opt));
    arx_opt = w_data(:,fix(phi_opt/pi*180)+51);
    rx_gain_sector(tt) = arx_opt' * H_NB_time_evo * atx_opt;
    SINR_sector(tt) = abs(rx_gain_sector(tt))^(2) / (noise_power);
    capacity_sector(tt) = log2(1+SINR_sector(tt)); % bps/Hz
    
    %% Evaluation of capacity of NB Tracking
%     arx_opt_NBtrack = exp(1j*(0:Nr-1)'*pi*sin(angle_est_rad(tt+1)));
    arx_opt_NBtrack = w_data(:,fix(angle_est_rad(tt+1)/pi*180)+51);
%     arx_opt_NBtrack = exp(1j*(0:Nr-1)'*pi*sin(mean(ray_AOA_azim(tt,:))));
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
    U_est = U_sorted(:,1:L)*sqrt(Nr);
    V_est = V_sorted(:,1:L)*sqrt(Nt);

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
plot(t_range(1:20:end),angle_rx_range(rx_opt_index(2:20:end))/pi*180,'o','linewidth',2);hold on
plot(t_range(1:20:end),angle_est_rad(2:20:end)/pi*180,'x','linewidth',2);hold on
plot(t_range(1:20:end),mean(ray_AOA_azim(1:20:end,:),2)/pi*180,'linewidth',2);hold on
legend('Neighbor Angle Refinement','NB Tracking','True AOA Mean')
grid on

% %%  capacity plotting
% figure;
% plot(t_range, SINR_sector,'--o','linewidth',2);hold on
% plot(t_range, SINR_NBtrack,'--x','linewidth',2);hold on
% plot(t_range, SINR_trueCSI,'linewidth',2);hold on
% legend('Neighbor Angle Refinement','NB Tracking','Perfect CSI')
% grid on
%%
% figure;
% plot(t_range, 10*log10(SINR_sector),'-.','linewidth',2);hold on
% plot(t_range, 10*log10(SINR_NBtrack),'-.','linewidth',2);hold on
% plot(t_range, 10*log10(SINR_trueCSI),'linewidth',2);hold on
% % plot(t_range, 10*log10(SINR_trueCSI),'linewidth',2);hold on
% 
% legend('Neighbor Angle Refinement','NB Tracking','Perfect CSI')
% title('Gain')
% xlabel('Time [s]')
% ylabel('Tx+Rx Gain [dB]')
% grid on
%%  Gain plotting

% figure;
% plot(t_range, 10*log10(SINR_sector./SINR_trueCSI),'linewidth',2);hold on
% plot(t_range, 10*log10(SINR_NBtrack./SINR_trueCSI),'linewidth',2);hold on
% % plot(t_range, 10*log10(SINR_trueCSI),'linewidth',2);hold on
% 
% legend('Neighbor Angle Refinement','NB Tracking')
% title('Gain')
% xlabel('Time [s]')
% ylabel('Gain Drop [dB]')
% grid on
%%
figure
[b,a] = ecdf(10*log10(SINR_sector));
plot(a,1-b,'linewidth',2);hold on
[b,a] = ecdf(10*log10(SINR_NBtrack));
plot(a,1-b,'linewidth',2);hold on
[b,a] = ecdf(10*log10(SINR_trueCSI));
plot(a,1-b,'linewidth',2);hold on
grid on
legend('Neighbor Angle Refinement','NB Tracking','True CSI')

%%
power_AOA_dB = 10*log10(power_AOA);
AOA_range = linspace(-90,90,256);
[X, Y] = meshgrid(t_range(1:20:end),AOA_range);
figure
surf(X,Y,power_AOA_dB(:,1:20:end)-max(max(power_AOA_dB(:,1:20:end))))
colorbar
colormap default
az = 0;
el = 90;
view(az, el);
% xlim([-90,30])
ylim([-10,20])
caxis([-20 0])
% title('Angular Power Profile (dB)')
xlabel('Time [s]')
ylabel('AoA [deg]')
