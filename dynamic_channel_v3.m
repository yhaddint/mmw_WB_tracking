% test whether it is possible to get reasonable "initial guess" of channle
% parameter
clear;clc;close all
rng(1)
%% other parameter
plot_ellipse = 1;
%% Goemetry Parameters
Nt = 64;
Nr = 32;
D_bs2ue = 40;
loc0_ue = [-D_bs2ue/2,0];
loc0_bs = [D_bs2ue/2,0];
speed_c = 3e8;
fc = 28e9;
%% Large Scale Parameter (Multiple Clusters)
cluster_delay = [300e-9,200e-9];
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
    raygain(2,:) = zeros(1,ray_num);
    
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

figure(99)
plot(loc0_ue(1),loc0_ue(2),'x');hold on
plot(loc0_bs(1),loc0_bs(2),'x');hold on
plot(loc_cluster_total(:,1),loc_cluster_total(:,2),'o');hold on
title('Geometric Stochastic Multipath')
axis([-80,80,-80,80])
xlabel('X Coordn (m)')
ylabel('Y Coordn (m)')
grid on



%% Frequency Domain MIMO Channel & Frequency Domain Equivalent Beamspace Channel
Nfft = 512;
H_freq0 = get_H_freq2( raygain(1,:),...
                       raydelay(1,:),...
                       rayAOA(1,:),...
                       rayAOD(1,:),...
                       1,...
                       ray_num,...
                       Nt, Nr, fc);

%%
phi1 = mean(rayAOA(1,:),2);
phi_round1 = round(mean(rayAOA(1,:),2)/pi*2^6)/2^6*pi;
arx1 = exp(1j*(0:Nr-1)'*pi*sin(phi_round1));
% phi_round1_leftasst = (round(mean(rayAOA(1,:),2)/pi*2^6)-0.25)/2^6*pi;
% arx1_leftasst = exp(1j*(0:Nr-1)'*pi*sin(phi_round1_leftasst));
% phi_round1_rightasst = (round(mean(rayAOA(1,:),2)/pi*2^6)+0.25)/2^6*pi;
% arx1_rightasst = exp(1j*(0:Nr-1)'*pi*sin(phi_round1_rightasst));

theta1 = mean(rayAOD(1,:),2);
atx1 = exp(1j*(0:Nt-1)'*pi*sin(theta1));
% for kk=1:Nfft
%     H_BB1(kk,1) = arx1'*squeeze(H_freq0(:,:,kk))*atx1;
% end

phi2 = mean(rayAOA(2,:),2);
phi_round2 = round(mean(rayAOA(1,:),2)/pi*2^6)/2^6*pi;
arx2 = exp(1j*(0:Nr-1)'*pi*sin(phi_round2));
theta2 = mean(rayAOD(2,:),2)
atx2 = exp(1j*(0:Nt-1)'*pi*sin(theta2));
%% quick look at post-beamforming wideband channel 
for kk=1:Nfft
    H_BB1(kk,1) = arx1'*squeeze(H_freq0(:,:,kk))*atx1;
end
figure
plot(10*log10(abs(H_BB1(:,1))))
grid on
xlabel('Subcarrier Index')
ylabel('Channel Gain (dB)')

%%
H_freq = get_H_freq2(raygain, raydelay, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr, fc);
rho = 1;
speed_v = [0,0];
t_range = (0:4:40)*1e-3;
rayAOA0 = rayAOA;
rayAOA_ez = repmat(mean(rayAOA0,2),1,ray_num);
raydelay0 = raydelay;
raygain0 = raygain;
for tt = 1:length(t_range)
    if tt==1
        raydelay_prev = ones(1,20)*300*1e-9;
        rayAOA_prev = ones(1,20)*phi1;
        raygain_prev = raygain;
    elseif tt>1
        raydelay_prev = tau_est;
        rayAOA_prev = AOA_est;
        raygain_prev = raygain;
%         raygain0 = alpha_est;
    end

    for kk=1:Nfft
        H_BB1(kk,tt) = arx1'*squeeze(H_freq(:,:,kk))*atx1;
        H_BB1_bf1(kk,tt) = arx1'*squeeze(H_freq(:,:,kk))*atx1;
%         H_BB1_bf2(kk,tt) = arx1_leftasst'*squeeze(H_freq(:,:,kk))*atx1;
%         H_BB1_bf3(kk,tt) = arx1_rightasst'*squeeze(H_freq(:,:,kk))*atx1;
    end

    % Alternative estimation parameters
    if tt>0
        %---------------------------------------------
        % tau_bar estimation with previous alpha and phi
        %---------------------------------------------
        for kk=1:Nfft
            bigPhi = sin(rayAOA_prev(1,:))-sin(phi_round1);
            bigTheta = -sin(rayAOD(1,:))+sin(theta1);
            BigTau(kk,:) = (-1j*2*pi/(Nfft*1e-9)*kk)...
            .*exp(-1j*2*pi*kk*raydelay_prev(1,:)/(1e-9*Nfft))...
            .*raygain(1,:)...    
            .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
            .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta));
            
            BigTau_const(kk,:) = exp(-1j*2*pi*kk*raydelay_prev(1,:)/(1e-9*Nfft))...
            .*raygain(1,:)...    
            .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
            .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta));
        end
%         H_BB1_pred(:,tt) = sum(BigTau_const,2)+BigTau*((raydelay(1,:)-raydelay_prev(1,:)).');
        deltataubar_est = pinv([real(sum(BigTau,2));imag(sum(BigTau,2))])...
            *[real(H_BB1(:,tt)-sum(BigTau_const,2));imag(H_BB1(:,tt)-sum(BigTau_const,2))];
        for kk=1:Nfft
        H_BB1_pred(kk,tt) = sum(raygain(1,:)...
                .*exp(-1j*2*pi*kk*(raydelay_prev(1,:)+deltataubar_est)/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
                .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta)));
        end

        
        %---------------------------------------------
        % phi_bar estimation using previous alpha and tau
        %---------------------------------------------
        for kk=1:Nfft
            bigPhi = sin(rayAOA_prev(1,:))-sin(phi_round1);
            bigTheta = -sin(rayAOD(1,:))+sin(theta1);
            BigPhi_const(kk,:) = raygain(1,:)...
                .*exp(-1j*2*pi*kk*(raydelay_prev(1,:)+deltataubar_est)/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
                .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta));
            for ray_index=1:ray_num
                phi_bar = phi_round1;
                phi_chan = rayAOA_prev(1,ray_index);
                BigPHI(kk,ray_index) = raygain(1,ray_index)...
                .*exp(-1j*2*pi*kk*(raydelay_prev(1,ray_index)+deltataubar_est)/(1e-9*Nfft))...
                .*((-1j*pi*Nr*cos(rayAOA_prev(1,ray_index)))*exp(1j*pi*Nr*bigPhi(ray_index))*(1-exp(1j*pi*bigPhi(ray_index)))...
                -(1-exp(1j*pi*Nr*bigPhi(ray_index)))*exp(1j*pi*bigPhi(ray_index))*(-1j*pi*cos(rayAOA_prev(1,ray_index))))...
                ./((1-exp(1j*pi*bigPhi(ray_index)))^2)...
                .*(1-exp(1j*pi*Nt*bigTheta(ray_index)))./(1-exp(1j*pi*bigTheta(ray_index)));
            end
        end
%         H_BB1_pred(:,tt) = sum(BigPhi_const,2)+BigPHI*((rayAOA(1,:)-rayAOA_prev(1,:)).');
        deltaAOAbar_est = pinv([real(sum(BigPHI,2));imag(sum(BigPHI,2))])...
            *([real(H_BB1(:,tt)-sum(BigPhi_const,2)); imag(H_BB1(:,tt)-sum(BigPhi_const,2))]);
%         AOA_est = rayAOA_prev(1,:) + deltaAOA_est.';
        
        
        %---------------------------------------------
        % tau estimation with previous alpha and phi
        %---------------------------------------------
%         deltaAOAbar_est = 0;
%         rayAOA_prev(1,:) = rayAOA_prev(1,:)+deltaAOAbar_est;
        for kk=1:Nfft
            bigPhi = sin(rayAOA_prev(1,:))-sin(phi_round1);
            bigTheta = -sin(rayAOD(1,:))+sin(theta1);
            BigTau(kk,:) = (-1j*2*pi/(Nfft*1e-9)*kk)...
            .*exp(-1j*2*pi*kk*(raydelay_prev(1,:))/(1e-9*Nfft))...
            .*raygain(1,:)...    
            .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
            .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta));
            
            BigTau_const(kk,:) = exp(-1j*2*pi*kk*(raydelay_prev(1,:))...
            /(1e-9*Nfft))...
            .*raygain(1,:)...    
            .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
            .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta));
        end
%         H_BB1_pred(:,tt) = sum(BigTau_const,2)+BigTau*((raydelay(1,:)-raydelay_prev(1,:)).');
        deltatau_est = pinv([real(BigTau);imag(BigTau)])...
            *[real(H_BB1(:,tt)-sum(BigTau_const,2));imag(H_BB1(:,tt)-sum(BigTau_const,2))];
        tau_est = raydelay_prev(1,:) + deltatau_est.';
%         
        %---------------------------------------------
        % Alpha estimation using previous tau and phi
        %---------------------------------------------
%         for kk=1:Nfft
%             bigPhi = sin(rayAOA_prev(1,:))-sin(phi_round1);
%             bigTheta = -sin(rayAOD(1,:))+sin(theta1);
%             BigAlpha(kk,:) = exp(-1j*2*pi*kk*(raydelay_prev(1,:)+deltatau_est.')/(1e-9*Nfft))...
%                 .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
%             .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta));
%         end
% %     H_BB1_pred(:,tt) = BigAlpha*(raygain(1,:).');
%         alpha_est = (pinv(BigAlpha)*H_BB1(:,tt)).';
% %     

        %---------------------------------------------
        % phi estimation using previous alpha and tau
        %---------------------------------------------
%         rayAOA_prev(1,:) = rayAOA_prev(1,:) + deltaAOA_est;
        for kk=1:Nfft
            bigPhi = sin(rayAOA_prev(1,:)+deltaAOAbar_est)-sin(phi_round1);
            bigTheta = -sin(rayAOD(1,:))+sin(theta1);
            BigPhi_const(kk,:) = raygain(1,:)...
                .*exp(-1j*2*pi*kk*(raydelay_prev(1,:)+deltataubar_est)/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
                .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta));
            for ray_index=1:ray_num
                phi_bar = phi_round1;
                phi_chan = rayAOA_prev(1,ray_index);
                BigPHI(kk,ray_index) = raygain(1,ray_index)...
                .*exp(-1j*2*pi*kk*(raydelay_prev(1,ray_index)+deltataubar_est)/(1e-9*Nfft))...
                .*((-1j*pi*Nr*cos(rayAOA_prev(1,ray_index)))*exp(1j*pi*Nr*bigPhi(ray_index))*(1-exp(1j*pi*bigPhi(ray_index)))...
                -(1-exp(1j*pi*Nr*bigPhi(ray_index)))*exp(1j*pi*bigPhi(ray_index))*(-1j*pi*cos(rayAOA_prev(1,ray_index))))...
                ./((1-exp(1j*pi*bigPhi(ray_index)))^2)...
                .*(1-exp(1j*pi*Nt*bigTheta(ray_index)))./(1-exp(1j*pi*bigTheta(ray_index)));
            end
        end
%         H_BB1_pred(:,tt) = sum(BigPhi_const,2)+BigPHI*((rayAOA(1,:)-rayAOA_prev(1,:)).');
        deltaAOA_est = pinv([real(BigPHI);imag(BigPHI)])...
            *[real(H_BB1(:,tt)-sum(BigPhi_const,2));imag(H_BB1(:,tt)-sum(BigPhi_const,2))];
        AOA_est = rayAOA_prev(1,:) + deltaAOA_est.';
        
        %---------------------------------------------
        % phi estimation using two assisted receiving beams
        %---------------------------------------------
%         rayAOA_prev(1,:) = rayAOA_prev(1,:) + deltaAOA_est;
%         for kk=1:Nfft
%             bigPhi_bf1 = sin(rayAOA_prev(1,:))-sin(phi_round1);
% %             bigPhi_bf2 = sin(rayAOA_prev(1,:))-sin(phi_round1_leftasst);
% %             bigPhi_bf3 = sin(rayAOA_prev(1,:))-sin(phi_round1_rightasst);
%             
%             bigTheta = -sin(rayAOD(1,:))+sin(theta1);
%             BigPhi_const_bf1(kk,:) = raygain(1,:)...
%                 .*exp(-1j*2*pi*kk*(raydelay_prev(1,:)+deltatau_est.')/(1e-9*Nfft))...
%                 .*(1-exp(1j*pi*Nr*bigPhi_bf1))./(1-exp(1j*pi*bigPhi_bf1))...
%                 .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta));
%             
% %             BigPhi_const_bf2(kk,:) = raygain(1,:)...
% %                 .*exp(-1j*2*pi*kk*(raydelay_prev(1,:)+deltatau_est.')/(1e-9*Nfft))...
% %                 .*(1-exp(1j*pi*Nr*bigPhi_bf2))./(1-exp(1j*pi*bigPhi_bf2))...
% %                 .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta));
% %             
% %             BigPhi_const_bf3(kk,:) = raygain(1,:)...
% %                 .*exp(-1j*2*pi*kk*(raydelay_prev(1,:)+deltatau_est.')/(1e-9*Nfft))...
% %                 .*(1-exp(1j*pi*Nr*bigPhi_bf3))./(1-exp(1j*pi*bigPhi_bf3))...
% %                 .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta));
%             for ray_index=1:ray_num
%                 phi_bar = phi_round1;
%                 phi_chan = rayAOA_prev(1,ray_index);
%                 
%                 BigPHI_bf1(kk,ray_index) = raygain(1,ray_index)...
%                 .*exp(-1j*2*pi*kk*(raydelay_prev(1,ray_index)+deltatau_est(ray_index))/(1e-9*Nfft))...
%                 .*((-1j*pi*Nr*cos(rayAOA_prev(1,ray_index)))*exp(1j*pi*Nr*bigPhi_bf1(ray_index))*(1-exp(1j*pi*bigPhi(ray_index)))...
%                 -(1-exp(1j*pi*Nr*bigPhi_bf1(ray_index)))*exp(1j*pi*bigPhi_bf1(ray_index))*(-1j*pi*cos(rayAOA_prev(1,ray_index))))...
%                 ./((1-exp(1j*pi*bigPhi_bf1(ray_index)))^2)...
%                 .*(1-exp(1j*pi*Nt*bigTheta(ray_index)))./(1-exp(1j*pi*bigTheta(ray_index)));
%             
% %                 BigPHI_bf2(kk,ray_index) = raygain(1,ray_index)...
% %                 .*exp(-1j*2*pi*kk*(raydelay_prev(1,ray_index)+deltatau_est(ray_index))/(1e-9*Nfft))...
% %                 .*((-1j*pi*Nr*cos(rayAOA_prev(1,ray_index)))*exp(1j*pi*Nr*bigPhi_bf2(ray_index))*(1-exp(1j*pi*bigPhi(ray_index)))...
% %                 -(1-exp(1j*pi*Nr*bigPhi_bf2(ray_index)))*exp(1j*pi*bigPhi_bf2(ray_index))*(-1j*pi*cos(rayAOA_prev(1,ray_index))))...
% %                 ./((1-exp(1j*pi*bigPhi_bf2(ray_index)))^2)...
% %                 .*(1-exp(1j*pi*Nt*bigTheta(ray_index)))./(1-exp(1j*pi*bigTheta(ray_index)));
% %                 
% %                 BigPHI_bf3(kk,ray_index) = raygain(1,ray_index)...
% %                 .*exp(-1j*2*pi*kk*(raydelay_prev(1,ray_index)+deltatau_est(ray_index))/(1e-9*Nfft))...
% %                 .*((-1j*pi*Nr*cos(rayAOA_prev(1,ray_index)))*exp(1j*pi*Nr*bigPhi_bf3(ray_index))*(1-exp(1j*pi*bigPhi(ray_index)))...
% %                 -(1-exp(1j*pi*Nr*bigPhi_bf3(ray_index)))*exp(1j*pi*bigPhi_bf3(ray_index))*(-1j*pi*cos(rayAOA_prev(1,ray_index))))...
% %                 ./((1-exp(1j*pi*bigPhi_bf3(ray_index)))^2)...
% %                 .*(1-exp(1j*pi*Nt*bigTheta(ray_index)))./(1-exp(1j*pi*bigTheta(ray_index)));
% %             
%             end
%         end
% %         H_BB1_pred(:,tt) = sum(BigPhi_const,2)+BigPHI*((rayAOA(1,:)-rayAOA_prev(1,:)).');
%         y1_LS_comp = H_BB1_bf1(:,tt)-sum(BigPhi_const_bf1,2);
%         A1_LS_comp  = BigPHI_bf1;
%         
% %         y2_LS_comp = H_BB1_bf2(:,tt)-sum(BigPhi_const_bf2,2);
% %         A2_LS_comp  = BigPHI_bf2;
% %         
% %         y3_LS_comp = H_BB1_bf3(:,tt)-sum(BigPhi_const_bf3,2);
% %         A3_LS_comp  = BigPHI_bf3;
% % %         
% % %         y_LS_3beam = [real(y1_LS_comp);imag(y1_LS_comp);
% % %                       real(y2_LS_comp);imag(y2_LS_comp)];
% % %                       real(y3_LS_comp);imag(y3_LS_comp)];
% %                   
% % %         A_LS_3beam = [real(A1_LS_comp);imag(A1_LS_comp);
% % %                       real(A2_LS_comp);imag(A2_LS_comp)];
% % %                       real(A3_LS_comp);imag(A3_LS_comp)];
%         
%         y_LS = [real(y1_LS_comp);imag(y1_LS_comp)];
%         A_LS = [real(A1_LS_comp);imag(A1_LS_comp)];
%         
% %         deltaAOA_est_3beam = pinv(A_LS_3beam)*y_LS_3beam;
%         deltaAOA_est = pinv(A_LS)*y_LS;
%         AOA_est = rayAOA_prev(1,:) + deltaAOA_est.';
%         
        %---------------------------------------------
        % Post-Beamforming Channel Using Estimated Parameter
        %---------------------------------------------
        bigPhi_est = sin(AOA_est)-sin(phi_round1);
        for kk=1:Nfft
        H_BB1_pred(kk,tt) = sum(raygain(1,:)...
                .*exp(-1j*2*pi*kk*tau_est/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi_est))./(1-exp(1j*pi*bigPhi_est))...
                .*(1-exp(1j*pi*Nt*bigTheta))./(1-exp(1j*pi*bigTheta)));
        end
        
        %---------------------------------------------
        % Pre-Beamforming MIMO Channel Using Estimated Parameter
        %---------------------------------------------
        % Est. chan
        H_freq_pred = get_H_freq2(raygain(1,:),...
                                  tau_est(1:ray_num),...
                                  AOA_est(1:ray_num),...
                                  rayAOD(1,1:ray_num),...
                                  1,...
                                  ray_num,...
                                  Nt, Nr, fc);
        
        % True chan.
        H_freq_true = get_H_freq2(raygain(1,:),...
                                  raydelay(1,:),...
                                  rayAOA(1,:),...
                                  rayAOD(1,:),...
                                  1,...
                                  ray_num,...
                                  Nt, Nr, fc);
        
        
        for kk=1:Nfft
            chan_error_pred(kk) = norm(squeeze(H_freq_pred(:,:,kk))-squeeze(H_freq_true(:,:,kk)),'fro')^2;
            chan_error_oldchan(kk) = norm(squeeze(H_freq0(:,:,kk))-squeeze(H_freq_true(:,:,kk)),'fro')^2;
            chan_pow(kk) = norm(squeeze(H_freq_true(:,:,kk)),'fro')^2;
        end


    end
    if tt==1
        H_BB1_pred(:,tt)=H_BB1(:,tt);
    end
    
%     if tt==8
%     figure
%     plot(abs(H_BB1(:,tt)));hold on
%     plot(abs(H_BB1_pred(:,tt)));hold on
%     plot(abs(H_BB1(:,1)));hold on
%     legend('chan now','chan prediction','chan prev')
%     end
    
    % pre-beamforming NMSE
    if tt==1
        MSE_oldchan(tt) = 0;
        MSE_pred(tt) = 0;
    else
        MSE_oldchan(tt) = mean(chan_error_oldchan)/mean(chan_pow);
        MSE_pred(tt) = mean(chan_error_pred)/mean(chan_pow);
    end
    
    
    % post-beamforming NMSE
    MSE1(tt) = norm(H_BB1(:,tt)-H_BB1(:,1),'fro')/norm(H_BB1(:,tt),'fro')
    MSE1_tracking(tt) = norm(H_BB1(:,tt)-H_BB1_pred(:,tt),'fro')/norm(H_BB1(:,tt),'fro')
    
%     MSE2(tt) = norm(H_BB2(:,tt)-H_BB2(:,1),'fro')/norm(H_BB2(:,tt),'fro')
%     MSE2_tracking(tt) = norm(H_BB2(:,tt)-H_BB2_pred(:,tt),'fro')/norm(H_BB2(:,tt),'fro')
end
%%
% total_pow = 10*log10(mean(abs(H_BB1).^2,1));
% figure;
% plot(total_pow-total_pow(1));hold on
% for tt = 1:length(t_range)
%     bigPhi = sin(phi_round1) - sin(mean(rayAOA0(1,:),2) + rayAOAdelta(1,tt));
%     BP(tt) = 20*log10(abs((1-exp(1j*pi*Nr*bigPhi))/(1-exp(1j*pi*bigPhi))));
% end
% plot(BP-BP(1));hold on

%% post-beamforming NMSE
figure
semilogy(1e3*t_range, (MSE1),'-o');hold on
semilogy(1e3*t_range, (MSE1_tracking),'-o');hold on
grid on
legend('w/o tracking','w/ tracking')
xlabel('Time (ms)')
ylabel('MSE')
title('Post-beamforming Chan. MSE')

%% pre-beamforming NMSE
figure
semilogy(1e3*t_range, (MSE_oldchan),'-o');hold on
semilogy(1e3*t_range, (MSE_pred),'-o');hold on
grid on
legend('w/o tracking','w/ tracking')
xlabel('Time (ms)')
ylabel('MSE')
title('pre-beamforming Chan. MSE')