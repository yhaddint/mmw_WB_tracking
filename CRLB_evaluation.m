%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;
rng(2); %random seed
% load probe_BF
%-------------------------------------
% System Parameters
%-------------------------------------
ray_num = 4; % Num of rays in a cluster
Nr = 16; % Number of antenna in Rx
M = 64; % Length of training
MCtimes = 200; % Num of Monte Carlo Sim.
AOAspread2 = (5/180*pi)^2;
AOAspread = sqrt(AOAspread2);

% For loop for Monte Carlo Simulations
for MCindex = 1:MCtimes
    
    % Receiver beamformer: 1) quasi-omni beam from random steering mtx; 2)
    % directional beam from angle steering vector
    probe_BF = (randi(2,Nr,M)*2-3) + 1j * (randi(2,Nr,M)*2-3);
    W = probe_BF./norm(probe_BF,'fro')*sqrt(Nr*M);
%     W = probe_BF(:,1:4)./norm(probe_BF(:,1:4),'fro');
    
    % Evaluate CRLB when the second ray is phi apart
    dphi_num = 100;
    dphi_range = linspace(2,15,dphi_num);
    
    % Evaluate CRLB when the second ray is phi apart
    SNR_num = 100;
    SNR_range = linspace(-15,30,SNR_num);
    
    % For loop for SNR
    for ss = 1:SNR_num
        
        % SNR
        sigman2 = 10^(-SNR_range(ss)/10);
        
        % AoA of rays with disired seperation
        phi = zeros(ray_num,1);
        phi0 = 0/180*pi;
        phi = phi0 + laprnd(ray_num,1,0,AOAspread);
%         randompm = randi(2,ray_num-1,1)*2-3;
%         phi = phi0 - AOAspread;
        
        % Pre-compute some vectors/matrices in FIM
        for rayindex = 1:ray_num
            
            % Complex gain
%             g_cmplx(rayindex) = (randn + 1j*randn)/sqrt(ray_num);
            g_cmplx = exp(1j*rand*2*pi)/ray_num;
%             g_cmplx = 1/ray_num;

            % Spatial response and its derivative over phi
            arx(:,rayindex) = exp(1j * pi * (0:Nr-1)' * sin(phi(rayindex)));
            Darx(:,rayindex) = 1j * pi * (0:Nr-1)' * cos(phi(rayindex));
        end
        
        % Evaluation of FIM with multiple rays
        
        % Pre-compute some equations, vecx, vecalpha, vecbeta
        % vecx is derivative of phi_0 in f(r,eta)
        vecx = zeros(M,1);
        for rr=1:ray_num
            vecx = vecx + g_cmplx*W'*(Darx(:,rr).*arx(:,rr));
        end
        
        % vecalpha is derivative of alpha in f(r,eta)
        vecalpha = zeros(M,1);
        for rr=1:ray_num
            vecalpha = vecalpha + W'*arx(:,rr);
        end
        
        % vecbeta is derivative of beta in f(r,eta)
        vecbeta = 1j * vecalpha;
        
        % Evaluation of J_{0,0} part, a 3 by 3 matrix
        % Partial of phi_0 and phi_0              
        J_00(1,1) = real(vecx' * vecx);

        % Partial of phi_0 and alpha
        J_00(1,2) = real(vecx' * vecalpha);

        % Partial of phi_r1 and beta_r2
        J_00(1,3) = real(vecx'*vecbeta);

        % Partial of alpha and phi_0
        J_00(2,1) = real(vecalpha'*vecx);

        % Partial of alpha and alpha
        J_00(2,2) = real(vecalpha'*vecalpha);

        % Partial of alpha and beta
        J_00(2,3) = real(vecalpha'*vecbeta);

        % Partial of beta and phi_0
        J_00(3,1) = real(vecbeta'* vecx);

        % Partial of beta and alpha
        J_00(3,2) = real(vecbeta'*vecalpha);

        % Partial of beta_r1 and beta_r2
        J_00(3,3) = real(vecbeta'*vecbeta);
        
        % Evaluate J_{0,Deltaphi} matrix
        for rr=1:ray_num
            
            % partial phi_0 and phi_r
            vec1 = vecx;
            vec2 = g_cmplx*W'*diag(Darx(:,rr))*arx(:,rr);
            J_0d(1,rr) = real(vec1'*vec2);
            
            % partial alpha and phi_r
            vec1 = vecalpha;
            vec2 = g_cmplx*W'*diag(Darx(:,rr))*arx(:,rr);
            J_0d(2,rr) = real(vec1'*vec2);
            
            % partial beta and phi_r
            vec1 = vecbeta;
            vec2 = g_cmplx*W'*diag(Darx(:,rr))*arx(:,rr);
            J_0d(3,rr) = real(vec1'*vec2);
            
        end
        % J(Deltaph,0) is transpose matrix of J_{0,Deltaphi}
        J_d0 = J_0d.';
        
        % Evaluate J{Deltaphi,Deltaphi} matrix
        for r1 = 1:ray_num
            for r2 = 1:ray_num
                vec1 = g_cmplx*W'*diag(Darx(:,r1))*arx(:,r1);
                vec2 = g_cmplx*W'*diag(Darx(:,r2))*arx(:,r2);
                J_dd(r1,r2) = real(vec1'*vec2);
            end
        end
        
        % Evaluate J matrix for apriori distribution of dphi
        J_p = [zeros(3,3), zeros(ray_num,3).';zeros(ray_num,3),diag(ones(ray_num,1)*sqrt(2)/AOAspread)];
        
        % Evaluate FIM
        J = 2/sigman2*[J_00, J_0d; J_d0, J_dd] + J_p;
        
        % RMSE evaluation from CRLB perspective
        % CRLB of first ray when there are multiple rays
        temp = inv(J);
        CRLB_multiple(ss,MCindex) = sqrt(temp(1,1))*(1/pi*180);
        CRLB_rest1(ss,MCindex) = sqrt(temp(4,4))*(1/pi*180);
        CRLB_rest2(ss,MCindex) = sqrt(temp(5,5))*(1/pi*180);
        
        
        % Evaluation of FIM with single rays
        temp = inv(J(1:3,1:3));
        % RMSE evaluation from CRLB perspective
        CRLB_single(ss,MCindex) = sqrt(temp(1,1))*(1/pi*180);
%         CRLB_single(pp,MCindex) = sqrt(1/(J(1,1) - J(1,2:3)*inv(J(2:3,2:3))*J(2:3,1)))*(1/pi*180);

    end
end
%% Plot CRLB of angle est./tracking
figure
semilogy(SNR_range,mean(CRLB_single,2));hold on
semilogy(SNR_range,mean(CRLB_multiple,2));hold on
grid on
legend('Single Ray','Multple Rays')

figure
semilogy(SNR_range,mean(CRLB_rest1,2));hold on
semilogy(SNR_range,mean(CRLB_rest1,2));hold on
grid on
legend('CRLB Ray 1','CRLB Ray 2')
