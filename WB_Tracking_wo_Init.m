% performance with perfect initial parameter...
clear;clc;close all
% rng(2)
%% other parameter
plot_ellipse = 0;
noise_pow = 1;

%% Goemetry Parameters
Nt = 16;
Nr = 16;
D_bs2ue = 40;
loc0_ue = [-D_bs2ue/2,0];
loc0_bs = [D_bs2ue/2,0];
speed_c = 3e8;
fc = 28e9;
%% Large Scale Parameter (Multiple Clusters)
cluster_delay = [300e-9,250e-9];
cluster_num = 1;
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
    loc_cluster = repmat(loc_cluster_cnt,ray_num,1)+[randn(ray_num,1)*0,randn(ray_num,1)*0];
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
SNR = -15:5:30;
errorSD_noisy = zeros(length(SNR),1); %RMSE
for noise_index = 1:length(SNR);
    angle_est = zeros(1000, 1);
    for i = 1:1000 %1000 Monte Carlo simulations
        y = 0; %received signal
        true_rayAOA = zeros(1, cluster_num);
        true_rayAOD = zeros(1, cluster_num);
        for cluster_index = 1:cluster_num
            delay = cluster_delay(cluster_index);
            ellipse_len = delay*speed_c;
            a = ellipse_len/2;
            c = D_bs2ue/2;
            b = sqrt(a^2-c^2);
            if cluster_index==1
                AOA = pi/2;
            else
                AOA = -pi/2;
            end
            loc_cluster_cnt = [a*cos(AOA),b*sin(AOA)];
            temp = loc_cluster_cnt - loc0_ue;
            true_rayAOA(1,cluster_index) = angle(temp(1)+1j*temp(2));
            temp = loc_cluster_cnt - loc0_bs;
            true_rayAOD(1,cluster_index) = angle(-temp(1)-1j*temp(2));
            atx = exp(1j*(0:Nt-1)'*pi*sin(true_rayAOD(1,cluster_index)))/sqrt(Nt); %we assume transmitter is perfectly aligned
            arx = exp(1j*(0:Nr-1)'*pi*sin(true_rayAOA(1,cluster_index)))/sqrt(Nr); %we assume transmitter is perfectly aligned
            arx1 = exp(1j*(0:Nr-1)'*pi*sin(true_rayAOA(1,cluster_index))+2/180*pi)/sqrt(Nr);
            arx2 = exp(1j*(0:Nr-1)'*pi*sin(true_rayAOA(1,cluster_index))-2/180*pi)/sqrt(Nr);
            arx3 = exp(1j*(0:Nr-1)'*pi*sin(true_rayAOA(1,cluster_index))-3/180*pi)/sqrt(Nr);
            arx4 = exp(1j*(0:Nr-1)'*pi*sin(true_rayAOA(1,cluster_index))+1/180*pi)/sqrt(Nr);
            %if we have multiple clusters, must save information about true
            %rayAOA, true rayAOD and steering vectors (e.g. make new vectors
            %and save them into these vectors)
            F = repmat(atx,1,40);%[atx atx atx atx]; %we assume that receiver makes 4 measurements (for example)
            W = (randi(2,Nr,40)*2-3)+1j*(randi(2,Nr,40)*2-3);
            y = y + W'*H_freq0(:,:,1)*F*ones(40,1); %we take just one subcarrier
            signal_power = norm(y)^2;
            noise_power = signal_power / 10^(SNR(noise_index)/10);
            noise_vec = (randn(40,1)*sqrt(noise_power/80) + 1i*randn(40,1)*sqrt(noise_power/80));
            y = y + noise_vec;
        end
        angle_est(i,1) = ml_angle(y, true_rayAOA, true_rayAOD, F, W, cluster_num, Nt, Nr);
    end
    errorSD_noisy(noise_index,1) = sqrt(mean((angle_est - ones(1000,1)*(true_rayAOA(1,1)/pi*180)).^2));
    %angle_error_var = var(angle_est);
end

close all
figure
semilogy(SNR, errorSD_noisy, 'o-r', 'Linewidth', 2)
xlabel('SNR [dB]')
ylabel('RMSE of AOA estimation')
title('RMSE vs. SNR[dB]')