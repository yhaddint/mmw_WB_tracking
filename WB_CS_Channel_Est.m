% baseline estimation and throughput
% reference: "frequency domain compressive channel estimation for frequency
% selective hybrid massive mimo systems"


%% parameters (same as in dynamic_channel_temp_v1)
close all; clear all
Nfft = 64;
Ts=1e-9*512/Nfft;

Nt = 16;
Nr = 16;
Nrf = 1; % # RF chains
Ns = Nrf;%
snr = 0; % dB
% t_range = (1:10:200)*1e-3;
max_iter = 2;

P_lin = 10^(snr/10);
noise_pow = 1;
ray_num = 20;
cluster_num = 2;
use_true_values = 0; % for debugging, set to 0 for actucal implementation
perfect_csi = 0; % 1 for sanity check, set to 0 for algorithm

%% get channel
t_range_new = [0:2.5:200]*1e-3;
dt = t_range_new(4)-t_range_new(3);
speed_rotate =  -100/180*pi;
capacity = zeros(max_iter, length(t_range_new));


for iter = 1: max_iter
    [H_freq,  Delta_freq, meanAOA, meanAOD, loc0_bs, loc0_ue, loc_cluster_total] = dynamic_channel_temp_v1(Nt, Nr, Nfft, Ts, ray_num, cluster_num );
    
    norm_factor = sqrt(mean(mean(mean(abs(H_freq).^2))));
    if(norm_factor==0)
        display('norm factor is 0');
        pause;
    end
    H_freq = H_freq /norm_factor ;
    
%     for t=t_range
%         H_freq_future{t==t_range} =  H_freq_future{t==t_range}/norm_factor;
%     end
    
    M = 200; % number of frames used in channel estimation
    
    %% compute  y(m)[k]
    F_tr = zeros(Nt, M);
    W_tr = zeros(Nr, M);
    y = zeros(Nfft,M);
    s = zeros(Nfft,M); % training symbols
%     s_dict = P_lin/Ns* [-1-1j, -1+1j, 1- 1j, 1+1j]; % QPSK symbols
    s_dict = P_lin/Ns* [1,1,1,1]; % QPSK symbols
    
    for m=1:M
        rand_phase = (2*pi).*rand(Nt,1);
        F_tr(:,m) = exp(1j*rand_phase)/sqrt(Nt);
        rand_phase = (2*pi).*rand(Nr,1);
        W_tr(:,m) = exp(1j*rand_phase)/sqrt(Nr);
        
        
        
        
        for k=1:Nfft
            s(k,m) = s_dict(randi(4,1));
            sig_only(k,m) = W_tr(:,m)' *  H_freq(:,:,k) *  F_tr(:,m) * s(k,m);
            noise_only(k,m) =W_tr(:,m)' * sqrt(noise_pow/2)*(randn(Nr,1) + 1j*randn(Nr,1));
            
            y(k,m) = W_tr(:,m)' *  H_freq(:,:,k) *  F_tr(:,m) * s(k,m) +  W_tr(:,m)' * sqrt(noise_pow/2)*(randn(Nr,1) + 1j*randn(Nr,1));
            
        end
    end
    
    snr_true = 10*log10(mean(abs(sig_only).^2, 2) ./ mean(abs(noise_only).^2, 2));
    
    %% generate dictionary
    % G= 32;
    Gt = 128;
    Gr = 128;
    fc = 60e9;
    
    A_t_grid = zeros(Nt,Gt);
    theta_grid = linspace(-pi/2,pi/2,Gt);
    
    A_r_grid = zeros(Nr,Gr);
    phi_grid = linspace(-pi/2,pi/2,Gr);
    
    for g = 1:Gt
        A_t_grid(:,g) = exp(1j*(0:Nt-1)'*pi*sin(theta_grid(g)));
    end
    for g=1:Gr
        A_r_grid(:,g) = exp(1j*(0:Nr-1)'*pi*sin(theta_grid(g)));
    end
    
    Psi = kron(conj(A_t_grid), A_r_grid);
    
    
    
    %since Nrf = 1, vec(y^(m)[k])is y^(m)[k]); so I am skipping vec ing
    Phi = zeros(M, Nt*Nr);
%     Phi_k_m  = zeros(M, Nt*Nr, Nfft);
    for m=1:M
        Phi(m,:) = kron( transpose(F_tr(:,m)),  (W_tr(:,m))' );
%         for k=1:Nfft
%             Phi_k_m(m,:,k) = kron(s(k,m)* transpose(F_tr(:,m)), W_tr(:,m)');
%         end
    end
    
    
    
    
    
    %% SW-OMP algorithm inputs: z[k], Phi, Psi, epsilon
    epsilon = noise_pow;
    C_w = zeros(M*Nrf, M*Nrf);
    
    for m=1:M
        C_w((m-1)*Nrf+1: m*Nrf , (m-1)*Nrf+1: m*Nrf) = W_tr(:,m)' *W_tr(:,m);
    end
    
    D_w = chol(C_w, 'upper');
    
    
    % Gamma = Phi*Psi;
    
    
    
    % A_t_true = exp(1j*(0:Nt-1)'*pi*sin(meanAOD*pi/180));
    % A_r_true = exp(1j*(0:Nr-1)'*pi*sin(meanAOA*pi/180));
    %
    % Psi_true = kron(conj(A_t_true), A_r_true);
    
    
%     for k=1:Nfft
        %     Gamma_true_k{k} =( squeeze(Phi_k_m(:,:,k)))*Psi_true;
        Gamma = Phi* Psi;
%     end
    
    
    %     m=2
    %     k=10
    
    %     y(k,m)
    %     temp_H = H_freq(:,:,k) ;
    %     kron(s(k,m)* transpose(F_tr(:,m)), W_tr(:,m)')*temp_H(:)
    
    % (squeeze(Phi_k_m(m,:,k)))*Psi_true* Delta_freq(k)
    
    % Gamma_true_k{k}(m)* Delta_freq(k)
    
    
    
    
    r = zeros(M,Nfft);
    
    for k=1:Nfft
        r(:,k) = transpose(y(k,:));
    end
    mse = inf;
    
    c = zeros(Gt*Gr,k);
    
    big_T = []; %common support
    
    big_T_all = 1:Gt*Gr;
    big_T_bar = big_T_all;
    
    L_est = 0; % size of big_T
    L_est_max = 1; % max number of rays =20;
    zeta = zeros(L_est_max,Nfft);
    num_iter = 0;
    
    % [~,p_true]=min(sum(abs(Psi-repmat(Psi_true,1,size(Psi,2))).^2));
    
    % while(mse>epsilon)
    % while(num_iter<100 && mse >1e-2)
    %     num_iter = num_iter + 1;
    for num_iter = 1:10
        for k=1:Nfft
%             c(:,k) = Gamma_k{k}' * (inv(D_w))' * r(:,k);
            c(:,k) = Gamma' *  r(:,k);
        end
        [~,p_star] = max(sum(abs(c),2));
        %     p_star=p_true;
        big_T = [big_T, p_star];
        L_est = L_est + 1   ;
        
        mse_k = zeros(1,Nfft);
        
        if(use_true_values==1)
            for k=1:Nfft
                zeta(1:L_est,k) = inv(Gamma_true_k{k}' * inv(C_w) * Gamma_true_k{k}) * Gamma_true_k{k}' * inv(C_w)  * transpose(y(k,:));
                
                % update residue
                r(:,k) =  transpose(y(k,:)) - Gamma_true_k{k}* zeta(1:L_est,k);
                
                mse_k(k) = r(:,k)' * inv(C_w)*r(:,k);
            end
        else
            
            for k=1:Nfft
                zeta(1:L_est,k) = pinv(Gamma(:,big_T)' * inv(C_w) * Gamma(:,big_T)) * Gamma(:, big_T)' * inv(C_w)  * transpose(y(k,:));
                
                % update residue
                r(:,k) =  transpose(y(k,:)) - Gamma(:,big_T)* zeta(1:L_est,k);
                
                mse_k(k) = r(:,k)' * inv(C_w)*r(:,k);
            end
        end
        mse = 1/(Nfft*M*Nrf) * sum(abs(mse_k));
        
    end
    
    %% reconstruct channel
    if(use_true_values)
        A_t_est = A_t_true;
        A_r_est = A_r_true;
    else
        AOAs = zeros(1,length(big_T));
        AODs = zeros(1,length(big_T));
        A_t_est = zeros(Nt, length(big_T));
        A_r_est = zeros(Nr, length(big_T));
        for p_cur=1:length(big_T)
            p_star = big_T(p_cur);
            
            theta_idx = max(floor(p_star/Gr),1);
            theta_est = theta_grid(theta_idx)*180/pi;
            A_t_est(:,p_cur) = [exp(1j*(0:Nt-1)'*pi*sin(theta_est*pi/180))];
            
            AODs(p_cur) = [theta_est];
            
            
            phi_idx = rem(p_star, Gr);
            if(phi_idx~=0)
                phi_est = phi_grid(phi_idx)*180/pi;
            else
                phi_est = phi_grid(1)*180/pi;
            end
            A_r_est(:,p_cur) = [exp(1j*(0:Nr-1)'*pi*sin(phi_est*pi/180))];
            AOAs(p_cur) = [phi_est];
            
            
        end
        
    end
    
    
    nmse_num = 0;
    nmse_den = 0;
    nmse = 0;
    %     sum(abs(zeta).^2, 2)
    
    
    %% Compute SINR and capacity
    
    L = 2;
    SINR_k = zeros(L,Nfft);
    t_idx = 1;
    
    
    
    for t_future= t_range_new
        if(t_future>0)
            % generate channel
             speed_v = [10,0];
                
                
                
                loc_ue = loc0_ue+speed_v*t_future;
                [raydelay, rayAOA_coord, rayAOD_coord ] = get_multipath(loc0_bs, loc_ue, loc_cluster_total,...
                    cluster_num, ray_num );
                
                
                rayAOA_coord = rayAOA_coord + speed_rotate * t_future;
                
                rayAOA_array = rayAOA_coord;
                meanAOA = mean(rayAOA_array,2)*180/pi;
                rayAOD_array = rayAOD_coord - pi;
                meanAOD = mean(rayAOD_array,2)*180/pi;
                
                
                H_new = get_H_freq_Shailesh( raydelay, rayAOA_array, rayAOD_array, cluster_num, ray_num, Nt, Nr, Nfft, Ts);
                H_new = H_new/ norm_factor;
        end
        
        for k = 1:Nfft
            if(t_future==0)
                vec_H_k_est =   kr_Shailesh(conj(A_t_est), (A_r_est))* (zeta(1:L_est,k) ./zeta(1:L_est,1)) ;%* sqrt(abs(zeta(1,k)).^2));

                
                H_k_est = reshape(vec_H_k_est,Nr, Nt);
                if(perfect_csi)
                    H_k_est = squeeze(H_freq(:,:,k));
                end
                
                
                
                [U,S,V] = svd(H_k_est);
                % arrange in descending order of S
                [diag_S, temp_idx] = sort(diag(S), 'descend');
                S_sorted  = diag(diag_S);
                U_sorted = U(:,temp_idx);
                V_sorted = V(:,temp_idx);
                
                
                U_est{k} = U_sorted(:,1:L);
                V_est{k} = V_sorted(:,1:L);
                
                
                H_k = squeeze(H_freq(:,:,k));
            else
                H_k = squeeze(H_new(:,:,k));
                

            end

            gain = U_est{k}' * H_k *V_est{k};
            
            if L==2
                SINR_k(1, k) = abs(gain(1,1))^(2) / ((abs(gain(1,2)).^2) + noise_pow);
                SINR_k(2, k) = abs(gain(2,2))^(2) / ((abs(gain(2,1)).^2) + noise_pow);
                
%                 SINR_k(l, k) = abs(gain(l,l))^2 / (sum(abs(gain(l,:)).^2) -abs(gain(l,l))^2  + noise_pow);
%                 SINR_k(l, k) = abs(gain(l,l))^2 / (sum(abs(gain(l,:)).^2) -abs(gain(l,l))^2  + noise_pow);
            elseif L==1
                SINR_k(1, k) = abs(gain(1,1))^(2) / (noise_pow);
            end
            
            [U_true,S_true,V_true] = svd(H_k);
            % arrange in descending order of S
            [diag_S_true, temp_idx] = sort(diag(S_true), 'descend');
            S_true_sorted  = diag(diag_S_true);
            U_true_sorted = U_true(:,temp_idx);
            V_true_sorted = V_true(:,temp_idx);
            
        end
        if L==1
            capacity(iter, t_future==t_range_new) = mean(log2(1+SINR_k)); % bpz/Hz
        elseif L==2
            capacity(iter, t_future==t_range_new) = mean(sum(log2(1+SINR_k))); % bpz/Hz
        end
    end
    
    fprintf('MC realization %d \n',iter)
    if(iter==100)
        display('100 iter completed');
    end
end


%%
figure(101)
plot(t_range_new, mean(capacity,1),'linewidth',3)
title(['Speed = ' num2str(speed_v(1)) ' m/s' ])
ylabel('Capacity(bps/Hz)')
xlabel('time (s)');
grid on;

% saveas(gcf,['baseline_capacity_' num2str(speed_v(1)) 'ms_random_scatterer_loc_K' num2str(Nfft)], 'fig')
% a=1;

