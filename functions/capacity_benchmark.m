function [capacity_refine,...
          capacity_NBtrack,...
          capacity_trueCSI] = capacity_benchmark(noise_pow,...
                                                    t_range,...
                                                    speed_rotate,...
                                                    plot_para,...
                                                    MCindex,...
                                                    L,...
                                                    beacon)
%CAPACITY_BENCHMARK_V1 Summary of this function goes here
%   Detailed explanation goes here

plot_ellipse = 0;
capacity = zeros(length(t_range),1);
print_stats = 0;
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
cluster_num = 2;
ray_num = 10;

%% Cluster Generation (Strongest Two?)
for cluster_index = 1:cluster_num
    delay = cluster_delay(cluster_index);
    ellipse_len = delay*speed_c;
    a = ellipse_len/2;
    c = D_bs2ue/2;
    b = sqrt(a^2-c^2);
    
%     AOA = rand*pi-pi/2;
    if cluster_index==1
        AOA = -pi/2;
    else
        AOA = pi/2;
    end

    loc_cluster_cnt = [a*cos(AOA),b*sin(AOA)];
    loc_cluster = repmat(loc_cluster_cnt,ray_num,1)+[randn(ray_num,1)*2,randn(ray_num,1)*2];
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
    
%     fprintf('Cluster %d:\n',cluster_index);
%     fprintf('DS Mean:    %4.2f ns\n', mean(raydelay(cluster_index,:)*1e9));
%     fprintf('DS Std Dev: %4.2f ns\n', sqrt(var(raydelay(cluster_index,:)*1e9)));
%     fprintf('AOAS Mean:    %4.2f degree\n', mean(rayAOA(cluster_index,:)/pi*180));
%     fprintf('AOAS Std Dev: %4.2f degree\n', sqrt(var(rayAOA(cluster_index,:)/pi*180)));
%     fprintf('AODS Mean:    %4.2f degree\n', mean(rayAOD(cluster_index,:)/pi*180));
%     fprintf('AODS Std Dev: %4.2f degree\n', sqrt(var(rayAOD(cluster_index,:)/pi*180)));
%     
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
%% Narrow Band Compressive Tracking Beacon (D.Ramasamy et el)
% each element is from \pm 1 \pm 1j
[~, M_beacon] = size(beacon);
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

%%
phi1 = mean(rayAOA(1,:),2);
phi_round1 = round(mean(rayAOA(1,:),2)/pi*2^6)/2^6*pi;
arx1 = exp(1j*(0:Nr-1)'*pi*sin(phi_round1))/sqrt(Nr);

phi_round1_leftasst = (round(mean(rayAOA(1,:),2)/pi*2^6)-0.25)/2^6*pi;
arx1_leftasst = exp(1j*(0:Nr-1)'*pi*sin(phi_round1_leftasst))/sqrt(Nr);
phi_round1_rightasst = (round(mean(rayAOA(1,:),2)/pi*2^6)+0.25)/2^6*pi;
arx1_rightasst = exp(1j*(0:Nr-1)'*pi*sin(phi_round1_rightasst))/sqrt(Nr);

theta1 = mean(rayAOD(1,:),2);
atx1 = exp(1j*(0:Nt-1)'*pi*sin(theta1))/sqrt(Nt);

phi2 = mean(rayAOA(2,:),2);
phi_round2 = round(mean(rayAOA(2,:),2)/pi*2^6)/2^6*pi;
arx2 = exp(1j*(0:Nr-1)'*pi*sin(phi_round2))/sqrt(Nr);
theta2 = mean(rayAOD(2,:),2);
atx2 = exp(1j*(0:Nt-1)'*pi*sin(theta2))/sqrt(Nt);
theta_est = [theta1;theta2];
%% quick look at post-beamforming wideband channel 
if abs(noise_pow - 10)<0.001
    NMSE = -10;
elseif abs(noise_pow - 1)<0.001
    NMSE = -20;
elseif abs(noise_pow - 0.1)<0.001
    NMSE = -30;
else
    fprintf('snr error');
end
H_error_scaler = sqrt(10^(NMSE/10));

for kk=1:Nfft
    H_error = (randn(Nr,Nt)+1j*randn(Nr,Nt))/sqrt(2)*H_error_scaler;
    H_BB1(kk,1) = arx1'*(squeeze(H_freq0(:,:,kk))+H_error)*atx1;
    H_BB2(kk,1) = arx2'*(squeeze(H_freq0(:,:,kk))+H_error)*atx2;
end
%% Beam Sweep Based Training
angle_tx_range = linspace(-pi/2,pi/2,Nt*2);
angle_rx_range = linspace(-pi/2,pi/2,Nr*2);

rx_RSS_bin = zeros(Nfft,1);
tx_index_exclusion = [];
rx_index_exclusion = [];
for ll = 1:L
    rx_RSS = ones(length(angle_rx_range),length(angle_tx_range))*(-100);
    for beam_tx_index = 1:length(angle_tx_range)
        if sum(beam_tx_index == tx_index_exclusion)==0
        theta = angle_tx_range(beam_tx_index);
        atx = exp(1j*(0:Nt-1)'*pi*sin(theta))/sqrt(Nt);
        for beam_rx_index = 1:length(angle_rx_range)
            if sum(beam_rx_index == rx_index_exclusion)==0
            phi = angle_rx_range(beam_rx_index);
            arx = exp(1j*(0:Nr-1)'*pi*sin(phi))/sqrt(Nr);
            for kk=1:Nfft
                AWGN = (randn + 1j*randn)/sqrt(2);
                rx_RSS_bin(kk) = arx'*(squeeze(H_freq0(:,:,kk))) * atx + AWGN;
            end
            rx_RSS(beam_rx_index,beam_tx_index) = 10*log10(sum(abs(rx_RSS_bin).^2));
            end
        end
        end
    end
    [~,tx_max_index(ll)] = max(max(rx_RSS));
    [~,rx_max_index(ll)] = max(max(transpose(rx_RSS)));
    tx_index_exclusion = tx_max_index(ll)-5:tx_max_index(ll)+5;
    rx_index_exclusion = rx_max_index(ll)-5:rx_max_index(ll)+5;
    
end

tx_angle = angle_tx_range(tx_max_index)/pi*180;
rx_angle = angle_rx_range(rx_max_index)/pi*180;
%
% figure
% surf(angle_tx_range/pi*180,angle_rx_range/pi*180,rx_RSS)
% xlabel('AOD [deg]')
% ylabel('AOA [deg]')


%% Channel change in BB1 and BB2
rho = 0.999;
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
    [raydelay, rayAOA, rayAOD ] = get_multipath_v2(loc0_bs,...
                                                   loc_ue,...
                                                   loc_cluster_total,...
                                                   cluster_num,...
                                                   ray_num,...
                                                   print_stats);
                                               
    rayAOA = rayAOA + speed_rotate * (tt-1) * dt;
%     raygain = raygain.*exp(1j*rand(cluster_num,ray_num)*2*pi*sqrt(1-rho^2));
    raygain = rho * raygain + (randn(cluster_num,ray_num)+1j*randn(cluster_num,ray_num))...
                                /sqrt(2)*sqrt(1-rho^2);
                            
%     H_freq = get_H_freq2(raygain, raydelay, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr);
%     H_freq = H_freq / norm_factor ;

    raygain_true_BB1(tt,:)  = raygain(1,:);
    raygain_true_BB2(tt,:)  = raygain(2,:);
    
    raydelay_true_BB1(tt,:)  = raydelay(1,:);
    raydelay_true_BB2(tt,:)  = raydelay(2,:);
    
    rayAOA_true_BB1(tt,:)  = rayAOA(1,:);
    rayAOA_true_BB2(tt,:)  = rayAOA(2,:);
end
rx_opt_index = zeros(length(t_range)+1,L);
for ll = 1:L
    rx_opt_index(1,ll) = rx_max_index(ll);
end
%% simulating mobile channel over time
for tt = 1:length(t_range)
    
    raygain(1,:) = raygain_true_BB1(tt,:);
    raydelay(1,:) = raydelay_true_BB1(tt,:);
    rayAOA(1,:) = rayAOA_true_BB1(tt,:);
    
    raygain(2,:) = raygain_true_BB2(tt,:);
    raydelay(2,:) = raydelay_true_BB2(tt,:);
    rayAOA(2,:) = rayAOA_true_BB2(tt,:);

    H_freq_raw = get_H_freq2(raygain, raydelay, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr);
    H_freq = H_freq_raw / norm_factor ;
%     MIMO_noise = (randn(Nr,Nt,Nfft) + 1j*randn(Nr,Nt,Nfft)) * sqrt(noise_pow/2);
%     H_freq_noisy = H_freq + MIMO_noise;
    
    %% Neighbor Angle Refinement
    refine_index_range = [0,1,-1,2,-2];
    refine_num = length(refine_index_range);
    for ll = 1:L
        theta_opt = angle_tx_range(tx_max_index(ll));
        atx_opt = exp(1j*(0:Nt-1)'*pi*sin(theta_opt))/sqrt(Nt);
        for refine_index=1:refine_num
            rx_refine_index = rx_opt_index(tt,ll) + refine_index_range(refine_index);
            if rx_refine_index>0 && rx_refine_index<length(angle_rx_range)+1
%                 fprintf('index is %4.2f\n',rx_refine_index);
                phi_refine = angle_rx_range(rx_refine_index);
                arx_refine = exp(1j*(0:Nr-1)'*pi*sin(phi_refine))/sqrt(Nr);
                for kk=1:Nfft
                    AWGN = (randn + 1j*randn)/sqrt(2);
                    rx_RSS_refine_bin(kk) = arx_refine'*(squeeze(H_freq(:,:,kk)))*atx_opt + AWGN;
                end
                rx_refine_RSS(refine_index) = 10*log10(sum(abs(rx_RSS_refine_bin).^2));
            else
                rx_refine_RSS(refine_index) = -100;
            end
                
        end
        [~,opt_index] = max(rx_refine_RSS);
        rx_opt_index(tt+1,ll) = rx_opt_index(tt,ll) + refine_index_range(opt_index);
    end
    %% NB tracking by updating angle & gain of each cluster
    %---------------------------------------------
    % Alpha estimation using previous Phi
    %---------------------------------------------
    theta_opt = theta1;
    atx_opt = exp(1j*(0:Nt-1)'*pi*sin(theta_opt))/sqrt(Nt);
    if tt==1
        phi_hat = [phi1; phi2];
    else
        phi_hat = phi_new;
    end
    kk = 1;
    rx_measure = zeros(M_beacon,1);
    Hk = squeeze(H_freq(:,:,kk));
    for m_beacon = 1:M_beacon
        rx_measure(m_beacon) = beacon(:,m_beacon)' * Hk * atx_opt;
    end
    arx_hat = zeros(Nr,L);
    for ll=1:L
        arx_hat(:,ll) = exp(1j*(0:Nr-1)'*pi*sin(phi_hat(ll)))/sqrt(Nr);
    end
    H_cal = beacon' * arx_hat;
    alpha_est = pinv(H_cal) * rx_measure;
    
    % quick watch received beacon
%     figure;
%     plot(abs(H_cal * alpha_est));hold on
%     plot(abs(rx_measure))
    
    %---------------------------------------------
    % Phi estimation using previous Alpha
    %---------------------------------------------
    dPhi = zeros(Nr,L);
    for ll=1:L
        for nn=1:Nr
            dPhi(nn,ll) = exp(1j*pi*(nn-1)*sin(phi_hat(ll)))/sqrt(Nr)...
                *(1j*pi*(nn-1)*cos(phi_hat(ll)));
        end
        D_cal(:,ll) = (beacon' * dPhi(:,ll)) * alpha_est(ll);
    end
    
    yr = rx_measure - H_cal * alpha_est;
    Q_cal = [real(D_cal);imag(D_cal)];
    tr = [real(yr);imag(yr)];
    
    dphi = pinv(Q_cal) * tr;
    phi_new = phi_hat + dphi;
    phi_array(:,tt) = phi_new;
    %% Evaluation of capacity of Neighbor Angle Refinement
    atx_opt = zeros(Nt,L);
    arx_opt = zeros(Nr,L);
    for ll=1:L
        theta_opt = angle_tx_range(tx_max_index(ll));
        atx_opt(:,ll) = exp(1j*(0:Nt-1)'*pi*sin(theta_opt))/sqrt(Nt);
        phi_opt = angle_rx_range(rx_opt_index(tt+1,ll));
        arx_opt(:,ll) = exp(1j*(0:Nr-1)'*pi*sin(phi_opt))/sqrt(Nr);
    end
    
    for kk=1:Nfft
        H_k = squeeze(H_freq(:,:,kk));
        gain = arx_opt' * H_k * atx_opt;
        if L==1
            SINR_refine(1, kk) = abs(gain(1,1))^(2) / (noise_pow);
        else
            SINR_refine(1, kk) = abs(gain(1,1))^(2) / ((abs(gain(1,2)).^2) + noise_pow);
            SINR_refine(2, kk) = abs(gain(2,2))^(2) / ((abs(gain(2,1)).^2) + noise_pow);
        end
    end
    if L==2
        capacity_refine(tt) = mean(sum(log2(1+SINR_refine))); % bps/Hz
    elseif L==1
        capacity_refine(tt) = mean(log2(1+SINR_refine)); % bps/Hz
    end
    %% Evaluation of capacity of NB Tracking
    chan_hat_NB = zeros(Nr,Nt);
    for ll = 1:L
        arx_hat(:,ll) = exp(1j*(0:Nr-1)'*pi*sin(phi_hat(ll)))/sqrt(Nr);
        atx_hat(:,ll) = exp(1j*(0:Nt-1)'*pi*sin(theta_est(ll)))/sqrt(Nt);
        chan_hat_NB = chan_hat_NB + alpha_est(ll) * arx_hat(:,ll) * atx_hat(:,ll)';
    end
    [U,S,V] = svd(chan_hat_NB);
    
    for kk=1:Nfft
        H_k = squeeze(H_freq(:,:,kk));
        gain_NBtrack = U(:,1:2)' * H_k * V(:,1:2);
        if L==1
            SINR_NBtrack(1, kk) = abs(gain_NBtrack(1,1))^(2) / (noise_pow);
        else
            SINR_NBtrack(1, kk) = abs(gain_NBtrack(1,1))^(2) / ((abs(gain_NBtrack(1,2)).^2) + noise_pow);
            SINR_NBtrack(2, kk) = abs(gain_NBtrack(2,2))^(2) / ((abs(gain_NBtrack(2,1)).^2) + noise_pow);
        end
    end
    
    if L==2
        capacity_NBtrack(tt) = mean(sum(log2(1+SINR_NBtrack))); % bps/Hz
    elseif L==1
        capacity_NBtrack(tt) = mean(log2(1+SINR_NBtrack)); % bps/Hz
    end

    %% Benchmark with perfect CSI
    for kk=1:Nfft
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

    if plot_para
        figure

        plot(angle_rx_range(rx_opt_index(:,1))/pi*180,'linewidth',2);
        hold on
        plot(phi_array(1,:)/pi*180,'linewidth',2);hold on
        hold on
        plot(mean(rayAOA_true_BB1,2)/pi*180,'linewidth',2);
        hold on
        if L==2
            plot(angle_rx_range(rx_opt_index(:,2))/pi*180,'linewidth',2);
            hold on
            plot(phi_array(2,:)/pi*180,'linewidth',2);
            hold on
            plot(mean(rayAOA_true_BB2,2)/pi*180,'linewidth',2);
            hold on
            legend('Refinement1','NB Tracking1','True AOA1', 'Refinement2','NB Tracking2','True AOA2')
        else
            legend('Refinement1','NB Tracking1','True AOA1')
        end
        ylabel('AOA (deg)')
        titlename = ['MC',num2str(MCindex),', all AOA'];
        title(titlename)

    end

end

