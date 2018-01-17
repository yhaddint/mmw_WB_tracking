% 10/09/2017
% "Waterfall" figure of channel AOA & Gain. It shows similar pattern as
% work from Y.Wang et al (Huawei)

clear;clc;close all
% rng(2)
%% other parameter
plot_ellipse = 1;
noise_pow = 1;

%% Goemetry Parameters
Nt = 32;
Nr = 16;
D_bs2ue = 40;
loc0_ue = [-D_bs2ue/2,0];
loc0_bs = [D_bs2ue/2,0];
speed_c = 3e8;
fc = 28e9;
%% Large Scale Parameter (Multiple Clusters)
cluster_delay = [300e-9,250e-9];
cluster_num = 2;
ray_num = 20;

%% Cluster Generation (Strongest Two?)
for cluster_index = 1:cluster_num
    delay = cluster_delay(cluster_index);
    ellipse_len = delay*speed_c;
    a = ellipse_len/2;
    c = D_bs2ue/2;
    b = sqrt(a^2-c^2);
    
%     AOA = rand*pi-pi/2;
    if cluster_index==1
        AOA = pi/2;
    else
        AOA = -pi/2;
    end
    loc_cluster_cnt = [a*cos(AOA),b*sin(AOA)];
    loc_cluster = repmat(loc_cluster_cnt,ray_num,1)+[rand(ray_num,1)*4-2,rand(ray_num,1)*4-2];
    loc_cluster_total((cluster_index-1)*ray_num+1:cluster_index*ray_num,:) =  loc_cluster;
    
    for ray_index = 1:ray_num
        raydelay(cluster_index,ray_index) = (norm(loc_cluster(ray_index,:)-loc0_bs)...
                    +norm(loc_cluster(ray_index,:)-loc0_ue))/speed_c;
        temp = loc_cluster(ray_index,:) - loc0_ue;
        rayAOA(cluster_index,ray_index) = angle(temp(1)+1j*temp(2));
        temp = loc_cluster(ray_index,:) - loc0_bs;
        rayAOD(cluster_index,ray_index) = angle(-temp(1)-1j*temp(2));
    end
    raygain = exp(1j*rand(cluster_num, ray_num)*2*pi);
%     raygain(2,:) = zeros(1,ray_num);
    
    fprintf('Cluster %d:\n',cluster_index);
    fprintf('DS Mean:    %4.2f ns\n', mean(raydelay(cluster_index,:)*1e9));
    fprintf('DS Std Dev: %4.2f ns\n', sqrt(var(raydelay(cluster_index,:)*1e9)));
    fprintf('AOAS Mean:    %4.2f degree\n', mean(rayAOA(cluster_index,:)/pi*180));
    fprintf('AOAS Std Dev: %4.2f degree\n', sqrt(var(rayAOA(cluster_index,:)/pi*180)));
    fprintf('AODS Mean:    %4.2f degree\n', mean(rayAOD(cluster_index,:)/pi*180));
    fprintf('AODS Std Dev: %4.2f degree\n', sqrt(var(rayAOD(cluster_index,:)/pi*180)));
    
    % plot ellipse
    if plot_ellipse
        t = -pi:0.01:pi;
        x = a*cos(t);
        y = b*sin(t);
        figure(99)
        plot(x,y,'g--');hold on
    end
end

% plot of scattering ellipse
if plot_ellipse
    figure(99)
    plot(loc0_ue(1),loc0_ue(2),'x');hold on
    plot(loc0_bs(1),loc0_bs(2),'x');hold on
    plot(loc_cluster_total(:,1),loc_cluster_total(:,2),'o');hold on
    title('Geometric Stochastic Multipath')
    axis([-80,80,-80,80])
    xlabel('X Coordn (m)')
    ylabel('Y Coordn (m)')
    grid on
end
%% Frequency Domain MIMO Channel & Frequency Domain Equivalent Beamspace Channel
Nfft = 512;
H_freq0 = get_H_freq3( raygain,...
                       raydelay,...
                       rayAOA,...
                       rayAOD,...
                       cluster_num,...
                       ray_num,...
                       Nt, Nr);
norm_factor = sqrt(mean(mean(mean(abs(H_freq0).^2))));
% norm_factor = 1;
H_freq0 = H_freq0 / norm_factor ;

%%
phi1 = mean(rayAOA(1,:),2);
% phi_round1 = round(mean(rayAOA(1,:),2)/pi*2^6)/2^6*pi;
phi_round1 = phi1;
arx1 = exp(1j*(0:Nr-1)'*pi*sin(phi_round1))/sqrt(Nr);

phi_round1_leftasst = (round(mean(rayAOA(1,:),2)/pi*2^6)-0.25)/2^6*pi;
arx1_leftasst = exp(1j*(0:Nr-1)'*pi*sin(phi_round1_leftasst))/sqrt(Nr);
phi_round1_rightasst = (round(mean(rayAOA(1,:),2)/pi*2^6)+0.25)/2^6*pi;
arx1_rightasst = exp(1j*(0:Nr-1)'*pi*sin(phi_round1_rightasst))/sqrt(Nr);

theta1 = mean(rayAOD(1,:),2);
atx1 = exp(1j*(0:Nt-1)'*pi*sin(theta1))/sqrt(Nt);

phi2 = mean(rayAOA(2,:),2);
phi_round2 = round(mean(rayAOA(2,:),2)/pi*2^6)/2^6*pi;
phi_round2 = phi2;
arx2 = exp(1j*(0:Nr-1)'*pi*sin(phi_round2))/sqrt(Nr);
theta2 = mean(rayAOD(2,:),2);
atx2 = exp(1j*(0:Nt-1)'*pi*sin(theta2))/sqrt(Nt);


%% Channel change in BB1 and BB2
% rho = 1;
% speed_rotate = -200/180*pi;
% t_range = (4:4:200)*1e-3;
% power_delay = zeros(Nfft, length(t_range));
% for tt = 1:length(t_range)
%     t_now = t_range(tt);
%     dt = t_range(2)-t_range(1);
%     
%     % setup speed
%     if t_now<500e-3
%         speed_v = [5, 0];
%     else
%         speed_v = [5, 0];
%     end
%     
%     if tt==1
%         loc_ue = loc0_ue + speed_v * dt;
%     else
%         loc_ue = loc_ue + speed_v * dt;
%     end
%      
% %     clc
% %     fprintf('Time Evolution %2.4f s\n',t_range(tt));
%     [raydelay, rayAOA, rayAOD ] = get_multipath_v2(loc0_bs, loc_ue, loc_cluster_total,...
%                                             cluster_num, ray_num );
%     power_delay(:,tt) = analyze_delay_spead(raydelay, cluster_num, ray_num, Nfft);
%     
%     rayAOA = rayAOA + speed_rotate * (tt-1) * dt;
%     raygain = raygain.*exp(1j*rand(cluster_num,ray_num)*2*pi*sqrt(1-rho^2));
% 
%     raygain_true_BB1(tt,:)  = raygain(1,:);
%     raygain_true_BB2(tt,:)  = raygain(2,:);
%     
%     raydelay_true_BB1(tt,:)  = raydelay(1,:);
%     raydelay_true_BB2(tt,:)  = raydelay(2,:);
%     
%     rayAOA_true_BB1(tt,:)  = rayAOA(1,:);
%     rayAOA_true_BB2(tt,:)  = rayAOA(2,:);
% end
% 
% figure;
% surf(t_range,1:512,power_delay);
% hold on
% ylabel('Delay (ns)');
% xlabel('Simulation Time (s)');
% ylim([200,300])

%% 
% for tt = 1:length(t_range)
%     
%     raygain(1,:) = raygain_true_BB1(tt,:);
%     raydelay(1,:) = raydelay_true_BB1(tt,:);
%     rayAOA(1,:) = rayAOA_true_BB1(tt,:);
%     
%     raygain(2,:) = raygain_true_BB2(tt,:);
%     raydelay(2,:) = raydelay_true_BB2(tt,:);
%     rayAOA(2,:) = rayAOA_true_BB2(tt,:);
%     
%     H_freq = get_H_freq3(raygain, raydelay, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr);
%     H_freq = H_freq / norm_factor;  
% 
% end

%% Power-Delay and Power-AOA Profile Evolutions
speed_v = [-10,0];
t_range = (0:20:1000)*1e-3;
rho = 0.99;
for tt = 1:length(t_range)
    loc_ue = loc0_ue+speed_v*t_range(tt);

    [raydelay, rayAOA, rayAOD ] = get_multipath_v2(loc0_bs, loc_ue, loc_cluster_total,...
                                            cluster_num, ray_num );
    
    raygain = rho * raygain + (randn(cluster_num,ray_num)+1j*randn(cluster_num,ray_num))...
                                /sqrt(2)*sqrt(1-rho^2);
    
    clc
    fprintf('Time Evolution %2.4f s\n',t_range(tt));
    for cluster_index = 1 : cluster_num
    fprintf('Cluster %d:\n',cluster_index);
    fprintf('DS Mean:    %4.2f ns\n', mean(raydelay(cluster_index,:)*1e9));
    fprintf('DS Std Dev: %4.2f ns\n', sqrt(var(raydelay(cluster_index,:)*1e9)));
    fprintf('AOAS Mean:    %4.2f degree\n', mean(rayAOA(cluster_index,:)/pi*180));
    fprintf('AOAS Std Dev: %4.2f degree\n', sqrt(var(rayAOA(cluster_index,:)/pi*180)));
    fprintf('AODS Mean:    %4.2f degree\n', mean(rayAOD(cluster_index,:)/pi*180));
    fprintf('AODS Std Dev: %4.2f degree\n', sqrt(var(rayAOD(cluster_index,:)/pi*180)));
    end
    
    t_resolution = 4e-9;
    power_delay_raw(tt,:) = analyze_delay_spead(raygain, raydelay, cluster_num, ray_num, Nfft, t_resolution);
    H_freq = get_H_freq3(raygain, raydelay, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr);
    power_AOA(tt,:) = analyze_AOA_spead(H_freq, Nt, Nr, Nfft);
end
%%
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



