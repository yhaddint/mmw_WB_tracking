% clear;clc;
function [H_freq,  Delta_freq, meanAOA, meanAOD, loc0_bs, loc0_ue, loc_cluster_total] = dynamic_channel_temp_v1(Nt, Nr, Nfft, Ts, ray_num,cluster_num )
% rng('default')
%% other parameter
plot_ellipse = 1;
%% Goemetry Parameters
% Nt = 32;
% Nr = 32;
D_bs2ue = 60;
loc0_ue = [-D_bs2ue/2,0];
loc0_bs = [D_bs2ue/2,0];
speed_c = 3e8;
%Nfft = 1;

%% Large Scale Parameter (Multiple Clusters)
cluster_delay = [300e-9, 250e-9];

% 

%% Cluster Generation (Strongest Two?)
% AOA_ellipse_cluster_range_1 = linspace(pi/2-pi/6,  pi/2+pi/6, 10);
% AOA_ellipse_cluster_range_2 = linspace(-pi/2-pi/6,  -pi/2+pi/6, 10);
% AOA_ellipse_cluster = [datasample(AOA_ellipse_cluster_range_1,1), datasample(AOA_ellipse_cluster_range_2,1)];

for cluster_index = 1:cluster_num
    delay = cluster_delay(cluster_index);
    ellipse_len = delay*speed_c;
    a = ellipse_len/2;
    c = D_bs2ue/2;
    b = sqrt(a^2-c^2);
    
    %     AOA_ellipse = AOA_ellipse_cluster(cluster_index)
    if cluster_index==1
        AOA_ellipse = -pi/2;
    else
        AOA_ellipse = pi/2;
    end
%     AOA = rand*pi-pi/2;
    loc_cluster_cnt = [a*cos(AOA_ellipse),b*sin(AOA_ellipse)];
    loc_cluster = repmat(loc_cluster_cnt, ray_num,1)+[randn(ray_num,1)*2,randn(ray_num,1)*2];
    loc_cluster_total((cluster_index-1)*ray_num+1:cluster_index*ray_num,:) =  loc_cluster;
    
    
%     for ray_index = 1:ray_num
%         raydelay(cluster_index,ray_index) = (norm(loc_cluster(ray_index,:)-loc0_bs)...
%                     +norm(loc_cluster(ray_index,:)-loc0_ue))/speed_c;
%         temp = loc_cluster(ray_index,:) - loc0_ue;
%         rayAOA_coord(cluster_index,ray_index) = angle(temp(1)+1j*temp(2));
%         temp = loc_cluster(ray_index,:) - loc0_bs;
%         rayAOD_coord(cluster_index,ray_index) = angle(temp(1)+1j*temp(2));
%     end
    
%     fprintf('Cluster %d:\n',cluster_index);
%     fprintf('DS Mean:    %4.2f ns\n', mean(raydelay(cluster_index,:)*1e9));
%     fprintf('DS Std Dev: %4.2f ns\n', sqrt(var(raydelay(cluster_index,:)*1e9)));
%     fprintf('AOAS Mean:    %4.2f degree\n', mean(rayAOA_coord(cluster_index,:)/pi*180));
%     fprintf('AOAS Std Dev: %4.2f degree\n', sqrt(var(rayAOA_coord(cluster_index,:)/pi*180)));
%     fprintf('AODS Mean:    %4.2f degree\n', mean(rayAOD_coord(cluster_index,:)/pi*180));
%     fprintf('AODS Std Dev: %4.2f degree\n', sqrt(var(rayAOD_coord(cluster_index,:)/pi*180)));
    
    % plot ellipse
    if plot_ellipse
        t = -pi:0.01:pi;
        x = a*cos(t);
        y = b*sin(t);
        figure(99)
        plot(x,y,'g--');hold on
    end
end

    [raydelay, rayAOA_coord, rayAOD_coord]=  get_multipath(loc0_bs, loc0_ue, loc_cluster_total,...
        cluster_num, ray_num );

    
figure(99)
plot(loc0_ue(1),loc0_ue(2),'o');hold on
plot(loc0_bs(1),loc0_bs(2),'x');hold on
plot(loc_cluster_total(:,1),loc_cluster_total(:,2),'o');hold on
title('Geometric Stochastic Multipath')
axis([-80,80,-80,80])
grid on
% 
% figure
% subplot(311)
% hist(reshape(raydelay*1e9,cluster_num*ray_num,1),20)
% xlim([0,500])
% xlabel('Delay (ns)')
% title('Histogram of Multipath Delay')
% grid on
% subplot(312)
% hist(reshape(rayAOA/pi*180,cluster_num*ray_num,1),20)
% xlabel('AOA (degree)')
% xlim([-180,180])
% title('Histogram of Multipath AOA')
% grid on
% subplot(313)
% hist(reshape(rayAOD/pi*180,cluster_num*ray_num,1),20)
% xlabel('AOD (degree)')
% xlim([-180,180])
% title('Histogram of Multipath AOD')
% grid on

%% Frequency Domain MIMO Channel & Frequency Domain Equivalent Beamspace Channel
% meanAOA_coord = mean(rayAOA_coord,2)*180/pi ;

% if(mean(rayAOD_coord,2)<0)
%     meanAOD_array = mean(rayAOD_coord,2)*180/pi+180;
% else
%     meanAOD_array = 180-mean(rayAOD_coord,2)*180/pi;
% end


rayAOA_array = rayAOA_coord;
meanAOA = mean(rayAOA_array,2)*180/pi;
rayAOD_array = rayAOD_coord - pi;
meanAOD = mean(rayAOD_array,2)*180/pi;

[H_freq, Delta_freq] = get_H_freq_Shailesh( raydelay, rayAOA_array, rayAOD_array, cluster_num, ray_num, Nt, Nr, Nfft, Ts);
% [H_freq, Delta_freq] = get_H_freq2(raygain,  raydelay, rayAOA_array, rayAOD_array, cluster_num, ray_num, Nt, Nr, Nfft, Ts);

% phi = 120.17/180*pi;
% arx = exp(1j*(0:Nr-1)'*pi*sin(phi));
% theta = 168.68/180*pi;
% atx = exp(1j*(0:Nt-1)'*pi*sin(theta));
% for kk=1:Nfft
%     H_BB1(kk,1) = arx'*squeeze(H_freq(:,:,kk))*atx;
% end
% 
% % figure; plot(abs(H_BB1))
% %%

% speed_v = [16.67,0];
% t_iter = 1;
% 
% for t_future = t_range;
%     loc_ue = loc0_ue+speed_v*t_future;
%     [raydelay, rayAOA_coord, rayAOD_coord ] = get_multipath(loc0_bs, loc_ue, loc_cluster_total,...
%         cluster_num, ray_num );
%     
%     rayAOA_array = rayAOA_coord;
%     meanAOA = mean(rayAOA_array,2)*180/pi;
%     rayAOD_array = rayAOD_coord - pi;
%     meanAOD = mean(rayAOD_array,2)*180/pi;
% 
% 
%     H_freq_future{t_iter} = get_H_freq( raydelay, rayAOA_array, rayAOD_array, cluster_num, ray_num, Nt, Nr, Nfft, Ts);
%     t_iter = t_iter + 1;
% end
end
% 
% phi = 140.25/180*pi;
% arx = exp(1j*(0:Nr-1)'*pi*sin(phi));
% theta = 168.68/180*pi;
% atx = exp(1j*(0:Nt-1)'*pi*sin(theta));
% for kk=1:Nfft
%     H_BB1(kk,2) = arx'*squeeze(H_freq(:,:,kk))*atx;
% end
% figure(77)
% plot(abs(H_BB1(:,1:2)));hold on
% xlabel('subcarrier index')
% grid on
% norm(H_BB1(:,1)-H_BB1(:,2),'fro')/norm(H_BB1(:,1),'fro')
% 
% %% Power-Delay and Power-AOA Profile Evolutions
% speed_v = [3,0];
% t_range = (0:50:1000)*1e-3;
% for tt = 1:length(t_range)
%     loc_ue = loc0_ue+speed_v*t_range(tt);
%     clc
%     fprintf('Time Evolution %2.4f s\n',t_range(tt));
%     [raydelay, rayAOA, rayAOD ] = get_multipath(loc0_bs, loc_ue, loc_cluster_total,...
%                                             cluster_num, ray_num );
%     power_delay(tt,:) = analyze_delay_spead(raydelay, cluster_num, ray_num);
%     H_freq = get_H_freq( raydelay, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr);
%     power_AOA(tt,:) = analyze_AOA_spead(H_freq, Nt, Nr, Nfft);
% end
% 
% %% Power-AOA Profile Plotting
% power_delay = power_delay+1e-1;
% power_delay_dB = 10*log10(power_delay);
% delay_range = 1:512;
% [X, Y] = meshgrid(delay_range,t_range);
% figure
% subplot(211)
% surf(X,Y,power_delay_dB)
% colorbar
% colormap default
% az = 0;
% el = 90;
% view(az, el);
% xlim([240,260])
% title('Power Delay Profile Temporal Evolution for Cluster 1 (dB)')
% xlabel('Delay (ns)')
% ylabel('Simulation Time (ns)')
% subplot(212)
% surf(X,Y,power_delay_dB)
% colorbar
% colormap default
% view(az, el);
% xlim([340,360])
% title('Power Delay Profile Temporal Evolution for Cluster 2 (dB)')
% xlabel('Delay (ns)')
% ylabel('Simulation Time (ns)')
% 
% % Power-AOA Profile Plotting
% power_AOA_dB = 10*log10(power_AOA);
% AOA_range = linspace(-180,180,64);
% [X, Y] = meshgrid(AOA_range,t_range);
% figure
% surf(X,Y,power_AOA_dB)
% colorbar
% colormap default
% az = 0;
% el = 90;
% view(az, el);
% xlim([-180,180])
% title('Power AOA Profile Temporal Evolution for Clusters (dB)')
% xlabel('AOA (deg)')
% ylabel('Simulation Time (ns)')
% 
% 
% 
