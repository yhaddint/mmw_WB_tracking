function [ capacity]...
    = capacity_eval_v2(t_range, speed_rotate, plot_para, MCindex)
%CAPACITY_EVAL Summary of this function goes here
%   Detailed explanation goes here
plot_ellipse = 0;
capacity = zeros(length(t_range),1);
%% Goemetry Parameters
Nt = 16;
Nr = 16;
D_bs2ue = 40;
loc0_ue = [-D_bs2ue/2,0];
loc0_bs = [D_bs2ue/2,0];
speed_c = 3e8;
fc = 28e9;
%% Large Scale Parameter (Multiple Clusters)
cluster_delay = [300e-9,200e-9];
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
        AOA = -pi/3;
    else
        AOA = pi/3;
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
        rayAOD(cluster_index,ray_index) = angle(temp(1)+1j*temp(2));
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
%% quick look at post-beamforming wideband channel 
for kk=1:Nfft
    H_BB1(kk,1) = arx1'*squeeze(H_freq0(:,:,kk))*atx1;
    H_BB2(kk,1) = arx2'*squeeze(H_freq0(:,:,kk))*atx2;
end
% figure
% plot(10*log10(abs(H_BB1(:,1))))
% grid on
% xlabel('Subcarrier Index')
% ylabel('Channel Gain (dB)')

%% Channel change in BB1 and BB2
rho = 1;
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
    [raydelay, rayAOA, rayAOD ] = get_multipath(loc0_bs, loc_ue, loc_cluster_total,...
                                            cluster_num, ray_num );
    rayAOA = rayAOA + speed_rotate * (tt-1) * dt;
    raygain = raygain.*exp(1j*rand(cluster_num,ray_num)*2*pi*sqrt(1-rho^2));
%     H_freq = get_H_freq2(raygain, raydelay, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr);
%     H_freq = H_freq / norm_factor ;

    raygain_true_BB1(tt,:)  = raygain(1,:);
    raygain_true_BB2(tt,:)  = raygain(2,:);
    
    raydelay_true_BB1(tt,:)  = raydelay(1,:);
    raydelay_true_BB2(tt,:)  = raydelay(2,:);
    
    rayAOA_true_BB1(tt,:)  = rayAOA(1,:);
    rayAOA_true_BB2(tt,:)  = rayAOA(2,:);
end

%% Parameter tracking using gradient descent

raydelay_est_BB1 = zeros(length(t_range),ray_num);
rayAOA_est_BB1 = zeros(length(t_range),ray_num);
raygain_est_BB1 = zeros(length(t_range),ray_num);

raydelay_est_BB2 = zeros(length(t_range),ray_num);
rayAOA_est_BB2 = zeros(length(t_range),ray_num);
raygain_est_BB2 = zeros(length(t_range),ray_num);

for tt = 1:length(t_range)
    
    raygain(1,:) = raygain_true_BB1(tt,:);
    raydelay(1,:) = raydelay_true_BB1(tt,:);
    rayAOA(1,:) = rayAOA_true_BB1(tt,:);
    
    raygain(2,:) = raygain_true_BB2(tt,:);
    raydelay(2,:) = raydelay_true_BB2(tt,:);
    rayAOA(2,:) = rayAOA_true_BB2(tt,:);

    H_freq = get_H_freq2(raygain, raydelay, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr);
    H_freq = H_freq / norm_factor ;
    
    arx1 = exp(1j*(0:Nr-1)'*pi*sin(phi_round1))/sqrt(Nr);
    arx2 = exp(1j*(0:Nr-1)'*pi*sin(phi_round2))/sqrt(Nr);

    for kk=1:Nfft
        H_BB1(kk,tt) = arx1'*squeeze(H_freq(:,:,kk))*atx1;   
        H_BB2(kk,tt) = arx2'*squeeze(H_freq(:,:,kk))*atx2;   
    end
    
    if tt==1
        raydelay_est_BB1(tt,:) = raydelay(1,:);
        rayAOA_est_BB1(tt,:) = rayAOA(1,:);
        raygain_est_BB1(tt,:) = raygain(1,:)/norm_factor;
        
        raydelay_est_BB2(tt,:) = raydelay(2,:);
        rayAOA_est_BB2(tt,:) = rayAOA(2,:);
        raygain_est_BB2(tt,:) = raygain(2,:)/norm_factor;
    else
        raydelay_est_BB1(tt,:) = tau_est_BB1;
        rayAOA_est_BB1(tt,:) = phi_est_BB1;
        raygain_est_BB1(tt,:) = alpha_est_BB1;
        
        raydelay_est_BB2(tt,:) = tau_est_BB2;
        rayAOA_est_BB2(tt,:) = phi_est_BB2;
        raygain_est_BB2(tt,:) = alpha_est_BB2;
    end
    
    % gain is not updated in this algorithm
    alpha_est_BB1 = raygain_est_BB1(tt,:);
    alpha_est_BB2 = raygain_est_BB2(tt,:);

    % LMS update for channel parameters in BB1
    [tau_est_BB1, deltaAOA_est_BB1,  phi_est_BB1] ...
        = LMS_tracking( alpha_est_BB1,...
                        raydelay_est_BB1(tt,:),...
                        rayAOA_est_BB1(tt,:),...
                        phi_round1,...
                        rayAOD(1,:),...
                        theta1,...
                        H_BB1(:,tt),...
                        Nt, Nr, Nfft, ray_num, stepsize);

    phi_round1 = phi_round1 + stepsize * deltaAOA_est_BB1.';

        % LMS update for channel parameters in BB2
    [tau_est_BB2, deltaAOA_est_BB2,  phi_est_BB2] ...
        = LMS_tracking( alpha_est_BB2,...
                        raydelay_est_BB2(tt,:),...
                        rayAOA_est_BB2(tt,:),...
                        phi_round2,...
                        rayAOD(2,:),...
                        theta2,...
                        H_BB2(:,tt),...
                        Nt, Nr, Nfft, ray_num, stepsize);

    phi_round2 = phi_round2 + stepsize * deltaAOA_est_BB2.';
        %---------------------------------------------
        % Post-Beamforming Channel Using Estimated Parameter
        %---------------------------------------------
        bigPhi_est = sin(phi_est_BB1)-sin(phi_round1);
        bigTheta = -sin(rayAOD(1,:))+sin(theta1);
        for kk=1:Nfft
        H_BB1_pred(kk,tt) = sum(alpha_est_BB1(1,:)...
                .*exp(-1j*2*pi*kk*tau_est_BB1/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi_est))./(1-exp(1j*pi*bigPhi_est))...
                .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta)))...
                ./(Nt*Nr);
        end
        
        bigPhi_est = sin(phi_est_BB2)-sin(phi_round2);
        bigTheta = -sin(rayAOD(2,:))+sin(theta2);
        for kk=1:Nfft
        H_BB2_pred(kk,tt) = sum(alpha_est_BB2(1,:)...
                .*exp(-1j*2*pi*kk*tau_est_BB2/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi_est))./(1-exp(1j*pi*bigPhi_est))...
                .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta)))...
                ./(Nt*Nr);
        end

    
    % post-BF evaluation
    MSE_tracking_BB1(tt) = norm(H_BB1(:,tt)-H_BB1_pred(:,tt),'fro')/norm(H_BB1(:,tt),'fro');
    MSE_oldest_BB1(tt) = norm(H_BB1(:,tt)-H_BB1(:,1),'fro')/norm(H_BB1(:,tt),'fro');
    
    MSE_tracking_BB2(tt) = norm(H_BB2(:,tt)-H_BB2_pred(:,tt),'fro')/norm(H_BB2(:,tt),'fro');
    MSE_oldest_BB2(tt) = norm(H_BB2(:,tt)-H_BB2(:,1),'fro')/norm(H_BB2(:,tt),'fro');
    
    % capacity evaluation
    L = cluster_num;
    SINR_k = zeros(L,Nfft);
    noise_pow = 1;

    H_freq_est = get_H_freq2([alpha_est_BB1; alpha_est_BB2],...
                             [tau_est_BB1;   tau_est_BB2],...
                             [phi_est_BB1;   phi_est_BB2],...
                              rayAOD,...
                              cluster_num,...
                              ray_num,...
                              Nt, Nr);
    
    for kk=1:Nfft
        [U,S,V] = svd(H_freq_est(:,:,kk));
        [diag_S, temp_idx] = sort(diag(S), 'descend');
        S_sorted  = diag(diag_S);
        U_sorted = U(:,temp_idx);
        V_sorted = V(:,temp_idx); 
        U_est{kk} = U_sorted(:,1:L);
        V_est{kk} = V_sorted(:,1:L);
        
        % true chan
        H_k = squeeze(H_freq(:,:,kk));
        
        gain = U_est{kk}' * H_k *V_est{kk};
        
        for ll=1:L
            SINR_k(ll, kk) = abs(gain(ll,ll))^(2)...
                        / (sum(abs(gain(ll,:)).^2) -abs(gain(ll,ll))^2 + noise_pow);
        end
        
        capacity(tt) = mean(sum(log2(1+SINR_k))); % bps/Hz
    end
end

    if plot_para
        figure
        subplot(221)
        plot(1e9*mean(raydelay_est_BB1,2));
        hold on
        plot(1e9*mean(raydelay_true_BB1,2));
        legend('est.','true')
        ylabel('delay (ns)')
        titlename = ['MC',num2str(MCindex),', delay BB1'];
        title(titlename)

        subplot(222)
        plot(1e9*mean(raydelay_est_BB2,2));
        hold on
        plot(1e9*mean(raydelay_true_BB2,2));
        legend('est.','true')
        ylabel('delay (ns)')
        titlename = ['MC',num2str(MCindex),', delay BB2'];
        title(titlename)

        subplot(223)
        plot(mean(rayAOA_est_BB1,2)/pi*180);
        hold on
        plot(mean(rayAOA_true_BB1,2)/pi*180);
        legend('est.','true')
        ylabel('AOA (deg)')
        titlename = ['MC',num2str(MCindex),', AOA BB1'];
        title(titlename)

        subplot(224)
        plot(mean(rayAOA_est_BB2,2)/pi*180);
        hold on
        plot(mean(rayAOA_true_BB2,2)/pi*180);
        legend('est.','true')
        ylabel('AOA (deg)')
        titlename = ['MC',num2str(MCindex),', AOA BB2'];
        title(titlename)

    end
end

