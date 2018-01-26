% test whether it is possible to get reasonable "initial guess" of channle
% parameter
clear;clc;close all
rng(1)
%% other parameter
plot_ellipse = 0;
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
        AOA = pi/2;
    else
        AOA = -pi/2;
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
%% OMP approach 
G = 1000;
Nfft = 512;
tau_range = linspace(150,350,G)*1e-9;
krange = transpose(1:Nfft);
for gg=1:G
    Psi(:,gg) = exp(-1j*(2*pi)*krange*tau_range(gg)/(Nfft*1e-9));
end
s_tau = zeros(Nfft,ray_num);
h_res = H_BB1(:,1);
for rr=1:ray_num
    [~, kr_bb1(rr)] = max(abs(Psi'*h_res));
    s_tau(:,rr) = Psi(:,kr_bb1(rr));
    h_res = (eye(Nfft) - s_tau(:,1:rr)*pinv(s_tau(:,1:rr)))*H_BB1(:,1);
end

s_tau = zeros(Nfft,ray_num);
h_res = H_BB2(:,1);
for rr=1:ray_num
    [~, kr_bb2(rr)] = max(abs(Psi'*h_res));
    s_tau(:,rr) = Psi(:,kr_bb2(rr));
    h_res = (eye(Nfft) - s_tau(:,1:rr)*pinv(s_tau(:,1:rr)))*H_BB2(:,1);
end
% [sort(tau_range(kr)).' sort(raydelay(1,:)).']*1e9
%% Parameter refinement using Alternative Search in BB1
t_num = 50;
for tt = 1:t_num
    if tt==1
        rayAOA_prev = ones(1,ray_num)*phi1;
        raydelay_prev = tau_range(kr_bb1);
    elseif tt>1
        raydelay_prev = tau_est;
        rayAOA_prev = AOA_est;
        raygain_prev = alpha_est;
    end

    % Alternative estimation parameters
    
        %---------------------------------------------
        % Alpha estimation using previous tau and phi
        %---------------------------------------------
        for kk=1:Nfft
            bigPhi = sin(rayAOA_prev(1,:))-sin(phi_round1);
            bigTheta = -sin(rayAOD(1,:))+sin(theta1);
            BigAlpha(kk,:) = exp(-1j*2*pi*kk*(raydelay_prev(1,:))/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
                .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta))...
                ./(Nt*Nr);
        end
        alpha_est = (pinv(BigAlpha)*H_BB1(:,1)).';
%         H_BB1_pred(:,tt) = BigAlpha*(alpha_est.');
%         figure
%         plot(abs(H_BB1_pred(:,tt)));hold on
%         plot(abs(H_BB1(:,1)));hold on
%         legend('pred','true')
        
        
        %---------------------------------------------
        % tau estimation with previous alpha and phi
        %---------------------------------------------
        for kk=1:Nfft
            bigPhi = sin(rayAOA_prev(1,:))-sin(phi_round1);
            bigTheta = -sin(rayAOD(1,:))+sin(theta1);
            BigTau(kk,:) = (-1j*2*pi/(Nfft*1e-9)*kk)...
            .*exp(-1j*2*pi*kk*(raydelay_prev(1,:))/(1e-9*Nfft))...
            .*alpha_est...    
            .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
            .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta))...
            ./(Nt*Nr);
            
            BigTau_const(kk,:) = exp(-1j*2*pi*kk*(raydelay_prev(1,:))...
            /(1e-9*Nfft))...
            .*alpha_est...    
            .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
            .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta))...
            ./(Nt*Nr);
        end
%         H_BB1_pred(:,tt) = sum(BigTau_const,2)+BigTau*((raydelay(1,:)-raydelay_prev(1,:)).');
        deltatau_est = pinv([real(BigTau);imag(BigTau)])...
            *[real(H_BB1(:,1)-sum(BigTau_const,2));imag(H_BB1(:,1)-sum(BigTau_const,2))];
        tau_est = raydelay_prev(1,:) + deltatau_est.';
%         H_BB1_pred_real(:,tt) = real(BigTau) * deltatau_est + real(sum(BigTau_const,2));
%         H_BB1_pred_imag(:,tt) = imag(BigTau) * deltatau_est + imag(sum(BigTau_const,2));
%         figure
%         plot(abs(H_BB1_pred_real(:,tt)+1j*H_BB1_pred_imag(:,tt)));hold on
%         plot(abs(H_BB1(:,1)));hold on
%         legend('pred','true')

        %---------------------------------------------
        % phi estimation using previous alpha and tau
        %---------------------------------------------
%         rayAOA_prev(1,:) = rayAOA_prev(1,:) + deltaAOA_est;
        for kk=1:Nfft
            bigPhi = sin(rayAOA_prev(1,:))-sin(phi_round1);
            bigTheta = -sin(rayAOD(1,:))+sin(theta1);
            BigPhi_const(kk,:) = alpha_est(1,:)...
                .*exp(-1j*2*pi*kk*(tau_est)/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
                .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta))...
                ./(Nt*Nr);
            for ray_index=1:ray_num
                phi_bar = phi_round1;
                phi_chan = rayAOA_prev(1,ray_index);
                BigPHI(kk,ray_index) = alpha_est(ray_index)...
                .*exp(-1j*2*pi*kk*(tau_est(ray_index))/(1e-9*Nfft))...
                .*((-1j*pi*Nr*cos(rayAOA_prev(1,ray_index)))*exp(1j*pi*Nr*bigPhi(ray_index))*(1-exp(1j*pi*bigPhi(ray_index)))...
                -(1-exp(1j*pi*Nr*bigPhi(ray_index)))*exp(1j*pi*bigPhi(ray_index))*(-1j*pi*cos(rayAOA_prev(1,ray_index))))...
                ./((1-exp(1j*pi*bigPhi(ray_index)))^2)...
                .*(1-exp(1j*pi*Nt*bigTheta(ray_index)))./(1-exp(1j*pi*bigTheta(ray_index)))...
                ./(Nt*Nr);
            end
        end
%         H_BB1_pred(:,tt) = sum(BigPhi_const,2)+BigPHI*((rayAOA(1,:)-rayAOA_prev(1,:)).');
        deltaAOA_est = pinv([real(BigPHI);imag(BigPHI)])...
            *[real(H_BB1(:,1)-sum(BigPhi_const,2));imag(H_BB1(:,1)-sum(BigPhi_const,2))];
        AOA_est = rayAOA_prev(1,:) + deltaAOA_est.';
        
%         H_BB1_pred_real(:,tt) = real(BigPHI) * deltaAOA_est + real(sum(BigPhi_const,2));
%         H_BB1_pred_imag(:,tt) = imag(BigPHI) * deltaAOA_est + imag(sum(BigPhi_const,2));
%         figure
%         plot(abs(H_BB1_pred_real(:,tt)+1j*H_BB1_pred_imag(:,tt)));hold on
%         plot(abs(H_BB1(:,1)));hold on
%         legend('pred','true')
        
        %---------------------------------------------
        % Post-Beamforming Channel Using Estimated Parameter
        %---------------------------------------------
        bigPhi_est = sin(AOA_est)-sin(phi_round1);
        for kk=1:Nfft
        H_BB1_pred(kk,tt) = sum(alpha_est(1,:)...
                .*exp(-1j*2*pi*kk*tau_est/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi_est))./(1-exp(1j*pi*bigPhi_est))...
                .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta)))...
                ./(Nt*Nr);
        end

    MSE_init1(tt) = norm(H_BB1(:,1)-H_BB1_pred(:,tt),'fro')/norm(H_BB1(:,1),'fro');
end
% wrap-up everything for tracking
tau_est1 = tau_est;
AOA_est1= AOA_est;
alpha_est1 = alpha_est;
%% 
figure
plot(MSE_init1);hold on;grid on
xlabel('Iteration')
title('MSE in Initilization')
ylabel('MSE in BB1')
%% Parameter refinement using Alternative Search in BB2
t_num = 50;
for tt = 1:t_num
    if tt==1
        rayAOA_prev = ones(1,ray_num)*phi2;
        raydelay_prev = tau_range(kr_bb2);
    elseif tt>1
        raydelay_prev = tau_est;
        rayAOA_prev = AOA_est;
        raygain_prev = alpha_est;
    end

    % Alternative estimation parameters
    
        %---------------------------------------------
        % Alpha estimation using previous tau and phi
        %---------------------------------------------
        for kk=1:Nfft
            bigPhi = sin(rayAOA_prev(1,:))-sin(phi_round2);
            bigTheta = -sin(rayAOD(2,:))+sin(theta2);
            BigAlpha(kk,:) = exp(-1j*2*pi*kk*(raydelay_prev(1,:))/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
            .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta))...
            ./(Nt*Nr);
        end
        alpha_est = (pinv(BigAlpha)*H_BB2(:,1)).';
        H_BB2_pred(:,tt) = BigAlpha*(alpha_est.');
%         figure
%         plot(abs(H_BB1_pred(:,tt)));hold on
%         plot(abs(H_BB1(:,tt)));hold on
%         legend('pred','true')
        
        
        %---------------------------------------------
        % tau estimation with previous alpha and phi
        %---------------------------------------------
        for kk=1:Nfft
            bigPhi = sin(rayAOA_prev(1,:))-sin(phi_round2);
            bigTheta = -sin(rayAOD(2,:))+sin(theta2);
            BigTau(kk,:) = (-1j*2*pi/(Nfft*1e-9)*kk)...
            .*exp(-1j*2*pi*kk*(raydelay_prev(1,:))/(1e-9*Nfft))...
            .*alpha_est...    
            .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
            .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta))...
            ./(Nt*Nr);
            
            BigTau_const(kk,:) = exp(-1j*2*pi*kk*(raydelay_prev(1,:))...
            /(1e-9*Nfft))...
            .*alpha_est...    
            .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
            .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta))...
            ./(Nt*Nr);
        end
%         H_BB1_pred(:,tt) = sum(BigTau_const,2)+BigTau*((raydelay(1,:)-raydelay_prev(1,:)).');
        deltatau_est = pinv([real(BigTau);imag(BigTau)])...
            *[real(H_BB2(:,1)-sum(BigTau_const,2));imag(H_BB2(:,1)-sum(BigTau_const,2))];
        tau_est = raydelay_prev(1,:) + deltatau_est.';
%         H_BB2_pred_real(:,tt) = real(BigTau) * deltatau_est + real(sum(BigTau_const,2));
%         H_BB2_pred_imag(:,tt) = imag(BigTau) * deltatau_est + imag(sum(BigTau_const,2));
%         figure
%         plot(abs(H_BB2_pred_real(:,tt)+1j*H_BB2_pred_imag(:,tt)));hold on
%         plot(abs(H_BB2(:,1)));hold on
%         legend('pred','true')

        %---------------------------------------------
        % phi estimation using previous alpha and tau
        %---------------------------------------------
%         rayAOA_prev(1,:) = rayAOA_prev(1,:) + deltaAOA_est;
        for kk=1:Nfft
            bigPhi = sin(rayAOA_prev(1,:))-sin(phi_round2);
            bigTheta = -sin(rayAOD(2,:))+sin(theta2);
            BigPhi_const(kk,:) = alpha_est(1,:)...
                .*exp(-1j*2*pi*kk*(tau_est)/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
                .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta))...
                ./(Nt*Nr);
            for ray_index=1:ray_num
                phi_bar = phi_round1;
                phi_chan = rayAOA_prev(1,ray_index);
                BigPHI(kk,ray_index) = alpha_est(ray_index)...
                .*exp(-1j*2*pi*kk*(tau_est(ray_index))/(1e-9*Nfft))...
                .*((-1j*pi*Nr*cos(rayAOA_prev(1,ray_index)))*exp(1j*pi*Nr*bigPhi(ray_index))*(1-exp(1j*pi*bigPhi(ray_index)))...
                -(1-exp(1j*pi*Nr*bigPhi(ray_index)))*exp(1j*pi*bigPhi(ray_index))*(-1j*pi*cos(rayAOA_prev(1,ray_index))))...
                ./((1-exp(1j*pi*bigPhi(ray_index)))^2)...
                .*(1-exp(1j*pi*Nt*bigTheta(ray_index)))./(1-exp(1j*pi*bigTheta(ray_index)))...
                ./(Nt*Nr);
            end
        end
%         H_BB1_pred(:,tt) = sum(BigPhi_const,2)+BigPHI*((rayAOA(1,:)-rayAOA_prev(1,:)).');
        deltaAOA_est = pinv([real(BigPHI);imag(BigPHI)])...
            *[real(H_BB2(:,1)-sum(BigPhi_const,2));imag(H_BB2(:,1)-sum(BigPhi_const,2))];
        AOA_est = rayAOA_prev(1,:) + deltaAOA_est.';
        
%         H_BB1_pred_real(:,tt) = real(BigPHI) * deltaAOA_est + real(sum(BigPhi_const,2));
%         H_BB1_pred_imag(:,tt) = imag(BigPHI) * deltaAOA_est + imag(sum(BigPhi_const,2));
%         figure
%         plot(abs(H_BB2_pred_real(:,tt)+1j*H_BB2_pred_imag(:,tt)));hold on
%         plot(abs(H_BB2(:,1)));hold on
%         legend('pred','true')
        
        %---------------------------------------------
        % Post-Beamforming Channel Using Estimated Parameter
        %---------------------------------------------
        bigPhi_est = sin(AOA_est)-sin(phi_round2);
        for kk=1:Nfft
        H_BB2_pred(kk,tt) = sum(alpha_est(1,:)...
                .*exp(-1j*2*pi*kk*tau_est/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi_est))./(1-exp(1j*pi*bigPhi_est))...
                .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta)))...
                ./(Nt*Nr);
        end
%         figure
%         plot(abs(H_BB2_pred(:,tt)));hold on
%         plot(abs(H_BB2(:,1)));hold on
%         legend('pred','true')


    MSE_init2(tt) = norm(H_BB2(:,1)-H_BB2_pred(:,tt),'fro')/norm(H_BB2(:,1),'fro');
end

% wrap-up everything for tracking
tau_est2 = tau_est;
AOA_est2 = AOA_est;
alpha_est2 = alpha_est;

%%
% figure
% plot(MSE_init2);hold on;grid on
% xlabel('Iteration')
% title('MSE in Initilization')
% ylabel('MSE in BB2')
%% Channel change in BB1 and BB2
rho = 1;
t_range = (10:10:1000)*1e-3;
for tt = 1:length(t_range)
    t_now = t_range(tt);
    dt = t_range(2)-t_range(1);
    
    % setup speed
    if t_now<500e-3
        speed_v = [5, 0];
    else
        speed_v = [5, 0];
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

%% Parameter tracking using Alternative Search BB1

raydelay_est = zeros(length(t_range),ray_num);
rayAOA_est = zeros(length(t_range),ray_num);
raygain_est = zeros(length(t_range),ray_num);

for tt = 1:length(t_range)
    
    raygain(1,:) = raygain_true_BB1(tt,:);
    raygain(2,:) = raygain_true_BB2(tt,:);
    raydelay(1,:) = raydelay_true_BB1(tt,:);
    raydelay(2,:) = raydelay_true_BB2(tt,:);
    rayAOA(1,:) = rayAOA_true_BB1(tt,:);
    rayAOA(2,:) = rayAOA_true_BB2(tt,:);

    H_freq = get_H_freq2(raygain, raydelay, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr);
    H_freq = H_freq / norm_factor ;
    
    arx1 = exp(1j*(0:Nr-1)'*pi*sin(phi_round1))/sqrt(Nr);
    
    for kk=1:Nfft
        H_BB1(kk,tt) = arx1'*squeeze(H_freq(:,:,kk))*atx1;   
    end
    
    
    if tt==1
        raydelay_est(tt,:) = tau_est1;
        rayAOA_est(tt,:) = AOA_est1;
        raygain_est(tt,:) = alpha_est1;
%         raydelay_est(tt,:) = raydelay(1,:);
%         rayAOA_est(tt,:) = rayAOA(1,:);
%         raygain_est(tt,:) = raygain(1,:)/norm_factor;
    else
        raydelay_est(tt,:) = tau_est;
        rayAOA_est(tt,:) = AOA_est;
        raygain_est(tt,:) = alpha_est;
    end
    
    % Alternative estimation parameters
        
        alpha_est = raygain_est(tt,:);
        %---------------------------------------------
        % tau estimation with previous alpha and phi
        %---------------------------------------------
        for kk=1:Nfft
            bigPhi = sin(rayAOA_est(tt,:))-sin(phi_round1);
            bigTheta = -sin(rayAOD(1,:))+sin(theta1);
            BigTau(kk,:) = (-1j*2*pi/(Nfft*1e-9)*kk)...
            .*exp(-1j*2*pi*kk*(raydelay_est(tt,:))/(1e-9*Nfft))...
            .*alpha_est...    
            .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
            .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta))...
            ./(Nt*Nr);
            
            BigTau_const(kk,:) = exp(-1j*2*pi*kk*(raydelay_est(tt,:))/(1e-9*Nfft))...
            .*alpha_est...    
            .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
            .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta))...
            ./(Nt*Nr);
        end
        deltatau_est = pinv([real(BigTau);imag(BigTau)])...
            *[real(H_BB1(:,tt)-sum(BigTau_const,2));imag(H_BB1(:,tt)-sum(BigTau_const,2))];
        tau_est = raydelay_est(tt,:) + deltatau_est.';
        
        
        %---------------------------------------------
        % phi estimation using two assisted receiving beams
        %---------------------------------------------
        for kk=1:Nfft
            bigPhi_bf1 = sin(rayAOA_est(tt,:))-sin(phi_round1);
            
            bigTheta = -sin(rayAOD(1,:))+sin(theta1);
            BigPhi_const_bf1(kk,:) = alpha_est(1,:)...
                .*exp(-1j*2*pi*kk*tau_est(1,:)/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi_bf1))./(1-exp(1j*pi*bigPhi_bf1))...
                .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta))...
                ./(Nt*Nr);
            
            
            for ray_index=1:ray_num
                BigPHI_bf1(kk,ray_index) = alpha_est(ray_index)...
                .*exp(-1j*2*pi*kk*tau_est(ray_index)/(1e-9*Nfft))...
                .*(((1-exp(1j*pi*bigPhi_bf1(ray_index)))*(-exp(1j*pi*Nr*bigPhi_bf1(ray_index)))*(1j*pi*Nr*cos(rayAOA_est(ray_index))))...
                  -((1-exp(1j*pi*Nr*bigPhi_bf1(ray_index)))*(-exp(1j*pi*bigPhi_bf1(ray_index)))*(1j*pi*cos(rayAOA_est(ray_index)))))...
                ./((1-exp(-1j*pi*bigPhi_bf1(ray_index))).^2)...
                .*(1-exp(1j*pi*Nt*bigTheta(ray_index)))./(1-exp(1j*pi*bigTheta(ray_index)))...
                ./(Nt*Nr);
            
            end
        end
        
        y1_LS_comp = H_BB1(:,tt)-sum(BigPhi_const_bf1,2);
        A1_LS_comp  = BigPHI_bf1;
        
        y_LS = [real(y1_LS_comp);imag(y1_LS_comp)];
        A_LS = [real(A1_LS_comp);imag(A1_LS_comp)];
        
        deltaAOA_est = sum(A_LS,2)'*y_LS;
        AOA_est = rayAOA_est(tt,:) + (1e-7)*deltaAOA_est.';

        
        %---------------------------------------------
        % Post-Beamforming Channel Using Estimated Parameter
        %---------------------------------------------
        bigPhi_est = sin(AOA_est)-sin(phi_round1);
        for kk=1:Nfft
        H_BB1_pred(kk,tt) = sum(alpha_est(1,:)...
                .*exp(-1j*2*pi*kk*tau_est/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi_est))./(1-exp(1j*pi*bigPhi_est))...
                .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta)))...
                ./(Nt*Nr);
        end
        phi_round1 = phi_round1+(1e-7)*deltaAOA_est.';

    
    % post-BF evaluation
    MSE_tracking_BB1(tt) = norm(H_BB1(:,tt)-H_BB1_pred(:,tt),'fro')/norm(H_BB1(:,tt),'fro');
    MSE_oldest_BB1(tt) = norm(H_BB1(:,tt)-H_BB1(:,1),'fro')/norm(H_BB1(:,tt),'fro');
    
end
raydelay_est_BB1 = raydelay_est;
rayAOA_est_BB1 = rayAOA_est;
raygain_est_BB1 = raygain_est;

%% Parameter tracking using Alternative Search BB2

raydelay_est = zeros(length(t_range),ray_num);
rayAOA_est = zeros(length(t_range),ray_num);
raygain_est = zeros(length(t_range),ray_num);

for tt = 1:length(t_range)
    
    raygain(1,:) = raygain_true_BB1(tt,:);
    raygain(2,:) = raygain_true_BB2(tt,:);
    raydelay(1,:) = raydelay_true_BB1(tt,:);
    raydelay(2,:) = raydelay_true_BB2(tt,:);
    rayAOA(1,:) = rayAOA_true_BB1(tt,:);
    rayAOA(2,:) = rayAOA_true_BB2(tt,:);

    H_freq = get_H_freq2(raygain, raydelay, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr);
    H_freq = H_freq / norm_factor ;
    
    arx2 = exp(1j*(0:Nr-1)'*pi*sin(phi_round2))/sqrt(Nr);
    
    for kk=1:Nfft
        H_BB2(kk,tt) = arx2'*squeeze(H_freq(:,:,kk))*atx2;
    end
    
    
    if tt==1
        raydelay_est(tt,:) = tau_est2;
        rayAOA_est(tt,:) = AOA_est2;
        raygain_est(tt,:) = alpha_est2;
%         raydelay_est(tt,:) = raydelay(2,:);
%         rayAOA_est(tt,:) = rayAOA(2,:);
%         raygain_est(tt,:) = raygain(2,:)/norm_factor;
    else
        raydelay_est(tt,:) = tau_est;
        rayAOA_est(tt,:) = AOA_est;
        raygain_est(tt,:) = alpha_est;
    end
    
    % Alternative estimation parameters
    
        alpha_est = raygain_est(tt,:);
        
        %---------------------------------------------
        % tau estimation with previous alpha and phi
        %---------------------------------------------
        for kk=1:Nfft
            bigPhi = sin(rayAOA_est(tt,:))-sin(phi_round2);
            bigTheta = -sin(rayAOD(2,:))+sin(theta2);
            BigTau(kk,:) = (-1j*2*pi/(Nfft*1e-9)*kk)...
            .*exp(-1j*2*pi*kk*(raydelay_est(tt,:))/(1e-9*Nfft))...
            .*alpha_est...    
            .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
            .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta))...
            ./(Nt*Nr);
            
            BigTau_const(kk,:) = exp(-1j*2*pi*kk*(raydelay_est(tt,:))/(1e-9*Nfft))...
            .*alpha_est...    
            .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
            .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta))...
            ./(Nt*Nr);
        end
        deltatau_est = pinv([real(BigTau);imag(BigTau)])...
            *[real(H_BB2(:,tt)-sum(BigTau_const,2));imag(H_BB2(:,tt)-sum(BigTau_const,2))];
        tau_est = raydelay_est(tt,:) + deltatau_est.';
        

        
        %---------------------------------------------
        % phi estimation using two assisted receiving beams
        %---------------------------------------------
        for kk=1:Nfft
            bigPhi_bf1 = sin(rayAOA_est(tt,:))-sin(phi_round2);
            
            bigTheta = -sin(rayAOD(2,:))+sin(theta2);
            BigPhi_const_bf1(kk,:) = alpha_est...
                .*exp(-1j*2*pi*kk*tau_est(1,:)/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi_bf1))./(1-exp(1j*pi*bigPhi_bf1))...
                .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta))...
                ./(Nt*Nr);
            
            
            for ray_index=1:ray_num
                BigPHI_bf1(kk,ray_index) = alpha_est(ray_index)...
                .*exp(-1j*2*pi*kk*tau_est(ray_index)/(1e-9*Nfft))...
                .*(((1-exp(1j*pi*bigPhi_bf1(ray_index)))*(-exp(1j*pi*Nr*bigPhi_bf1(ray_index)))*(1j*pi*Nr*cos(rayAOA_est(ray_index))))...
                  -((1-exp(1j*pi*Nr*bigPhi_bf1(ray_index)))*(-exp(1j*pi*bigPhi_bf1(ray_index)))*(1j*pi*cos(rayAOA_est(ray_index)))))...
                ./((1-exp(-1j*pi*bigPhi_bf1(ray_index))).^2)...
                .*(1-exp(1j*pi*Nt*bigTheta(ray_index)))./(1-exp(1j*pi*bigTheta(ray_index)))...
                ./(Nt*Nr);
            end
        end
%         H_BB1_pred(:,tt) = sum(BigPhi_const,2)+BigPHI*((rayAOA(1,:)-rayAOA_prev(1,:)).');
        y1_LS_comp = H_BB2(:,tt)-sum(BigPhi_const_bf1,2);
        A1_LS_comp  = BigPHI_bf1;

        
        y_LS = [real(y1_LS_comp);imag(y1_LS_comp)];
        A_LS = [real(A1_LS_comp);imag(A1_LS_comp)];
        
        deltaAOA_est = sum(A_LS,2)'*y_LS;
        AOA_est = rayAOA_est(tt,:) + (1e-7)*deltaAOA_est.';
        
        %---------------------------------------------
        % Post-Beamforming Channel Using Estimated Parameter
        %---------------------------------------------
        bigPhi_est = sin(AOA_est)-sin(phi_round2);
        for kk=1:Nfft
        H_BB2_pred(kk,tt) = sum(alpha_est(1,:)...
                .*exp(-1j*2*pi*kk*tau_est/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi_est))./(1-exp(1j*pi*bigPhi_est))...
                .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta)))...
                ./(Nt*Nr);
        end
        phi_round2 = phi_round2+(1e-7)*deltaAOA_est.';

    
    % post-BF evaluation
    MSE_tracking_BB2(tt) = norm(H_BB2(:,tt)-H_BB2_pred(:,tt),'fro')/norm(H_BB2(:,tt),'fro');
    MSE_oldest_BB2(tt) = norm(H_BB2(:,tt)-H_BB2(:,1),'fro')/norm(H_BB2(:,tt),'fro');
    
end
raydelay_est_BB2 = raydelay_est;
rayAOA_est_BB2 = rayAOA_est;
raygain_est_BB2 = raygain_est;
%% post-beamforming NMSE
figure
plot(1e3*t_range, MSE_tracking_BB1,'linewidth',3);hold on;
plot(1e3*t_range, MSE_oldest_BB1,'linewidth',3);hold on
plot(1e3*t_range, MSE_tracking_BB2,'linewidth',3);hold on;
plot(1e3*t_range, MSE_oldest_BB2,'linewidth',3);hold on
grid on
xlabel('time (ms)')
title('MSE in Tracking')
ylabel('MSE')
legend('BB1 w/ tracking','BB1 w/o tracking','BB2 w/ tracking','BB2 w/o tracking')


%% reconstruct MIMO channel
L = cluster_num;
SINR_k = zeros(L,Nfft);
noise_pow = 1;
for tt=1:length(t_range)
    H_freq_est = get_H_freq2([raygain_est_BB1(tt,:);  raygain_est_BB2(tt,:)],...
                             [raydelay_est_BB1(tt,:); raydelay_est_BB2(tt,:)],...
                             [rayAOA_est_BB1(tt,:);   rayAOA_est_BB2(tt,:)],...
                              rayAOD,...
                              cluster_num,...
                              ray_num,...
                              Nt, Nr);
                          
    H_freq_true = get_H_freq2([raygain_true_BB1(tt,:);  raygain_true_BB2(tt,:)],...
                             [raydelay_true_BB1(tt,:); raydelay_true_BB2(tt,:)],...
                             [rayAOA_true_BB1(tt,:);   rayAOA_true_BB2(tt,:)],...
                              rayAOD,...
                              cluster_num,...
                              ray_num,...
                              Nt, Nr);
    H_freq_true = H_freq_true / norm_factor;                     
    
    
    for kk=1:Nfft
        [U,S,V] = svd(H_freq_est(:,:,kk));
        [diag_S, temp_idx] = sort(diag(S), 'descend');
        S_sorted  = diag(diag_S);
        U_sorted = U(:,temp_idx);
        V_sorted = V(:,temp_idx); 
        U_est{kk} = U_sorted(:,1:L);
        V_est{kk} = V_sorted(:,1:L);
        
        % true chan
        H_k = squeeze(H_freq_true(:,:,kk));
        
        gain = U_est{kk}' * H_k *V_est{kk};
        
        for ll=1:L
            SINR_k(ll, kk) = abs(gain(ll,ll))^(2)...
                        / (sum(abs(gain(ll,:)).^2) -abs(gain(ll,ll))^2 + noise_pow);
        end
        
        capacity(tt) = mean(sum(log2(1+SINR_k))); % bpz/Hz
    end
end
%%  capacity plotting
figure; plot(t_range, capacity)
        
