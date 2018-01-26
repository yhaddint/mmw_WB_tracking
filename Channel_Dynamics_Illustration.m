% 10/09/2017
% "Waterfall" figure of channel AOA & Gain. It shows similar pattern as
% work from Y.Wang et al (Huawei)


%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;close all
rng(2); %random seed
plot_ellipse = 0; % plot geolocations of Tx/Rx/Scatterers
print_stat = 1; % print channel parameter summary

%-------------------------------------
% System Parameters
%-------------------------------------
Nt = 32; % Number of Tx antennas (ULA)
Nr = 16; % Number of Tx antennas (ULA)
noise_pow = 1;
fc = 28e9; % Carrier frequency
cluster_num = 2; % Number of multipath clusters
ray_num = 20; % Number of intra-cluster rays
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

%% Power-Delay and Power-AOA Profile Evolutions
speed_v = [-10,0];
t_range = (0:20:1000)*1e-3;
rho = 0.99;
for tt = 1:length(t_range)
    
    loc_ue = loc0_ue+speed_v*t_range(tt);
    
    clc;
    fprintf('Time Evolution %2.4f s\n',t_range(tt));
    [raydelay, rayAOA, rayAOD ] = get_multipath_v2(loc0_bs,...
                                                   loc_ue,...
                                                   loc_cluster_total,...
                                                   cluster_num,...
                                                   ray_num,...
                                                   print_stat );
    
    raygain = rho * raygain + (randn(cluster_num,ray_num)+1j*randn(cluster_num,ray_num))...
                                /sqrt(2)*sqrt(1-rho^2);
    
    t_resolution = 4e-9;
    power_delay_raw(tt,:) = analyze_delay_spead(raygain, raydelay, cluster_num, ray_num, Nfft, t_resolution);
    H_freq = get_H_freq3(raygain, raydelay, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr);
    power_AOA(tt,:) = analyze_AOA_spead(H_freq, Nt, Nr, Nfft);
end

%% Figure plotting

AWGN = randn(length(t_range),Nfft)+1j*randn(length(t_range),Nfft);
power_delay = abs(power_delay_raw + AWGN * (1e-3)).^2;
power_delay_dB = 10*log10(power_delay);
delay_range = 1:Nfft;
[X, Y] = meshgrid(delay_range,t_range);

% Power-AOA & Power-Delay Profile Plotting
figure
% subplot(211)
surf(X,Y,power_delay_dB-3);hold on
colorbar
colormap default
lim = caxis;
caxis([-20 0])

az = 0;
el = 90;
view(az, el);
xlim([230,280])
title('Simulated Power Delay Profile, Omni-Directional (dB)')
xlabel('Delay (nano-sec)')
ylabel('Simulated Time (sec)')

% subplot(212)
% surf(X,Y,power_delay_dB)
% colorbar
% colormap default
% view(az, el);
% xlim([280,310])
% title('Power Delay Profile Temporal Evolution for Cluster 2 (dB)')
% xlabel('Delay (ns)')
% ylabel('Simulation Time (ns)')

% Power-AOA Profile Plotting
power_AOA_dB = 10*log10(power_AOA);
AOA_range = linspace(-90,90,64);
[X, Y] = meshgrid(AOA_range,t_range);
figure
surf(X,Y,power_AOA_dB-80)
colorbar
colormap default
az = 0;
el = 90;
view(az, el);
xlim([-90,30])
caxis([-20 0])
title('Simulated AOA v.s. Time (dB)')
xlabel('Rx Array Steering Angle (deg)')
ylabel('Simulated Time (sec)')



