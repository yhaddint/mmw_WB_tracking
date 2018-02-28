%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;
rng(2); %random seed
% load probe_BF
%-------------------------------------
% System Parameters
%-------------------------------------
ray_num = 20; % Num of rays in a cluster
Nr = 16; % Number of antenna in Rx
Nt = 32;
M = 16; % Length of training
MCtimes = 50; % Num of Monte Carlo Sim.
AOAspread2 = (5/180*pi)^2;
AOAspread = sqrt(AOAspread2);
AODspread = AOAspread;
SNR_num = 50;
SNR_range = linspace(-15,30,SNR_num);

% For loop for Monte Carlo Simulations (realization of g, angle spread, and beamformer)
for MCindex = 1:MCtimes
    
    clc
    fprintf('iteration %d:\n',MCindex);
    
    % Receiver beamformer: 1) quasi-omni beam from random steering mtx; 2)
    % directional beam from angle steering vector
    probe_Rx_BF = (randi(2,Nr,M)*2-3) + 1j * (randi(2,Nr,M)*2-3);
    W = probe_Rx_BF./norm(probe_Rx_BF,'fro')*sqrt(Nr*M);
    
    probe_Tx_BF = (randi(2,Nt,M)*2-3) + 1j * (randi(2,Nt,M)*2-3);
    F = probe_Tx_BF./norm(probe_Tx_BF,'fro')*sqrt(Nt*M);   

    % AoA of rays with disired seperation
    phi = zeros(ray_num,1);
    phi0 = 0/180*pi;
    phi = phi0 + laprnd(ray_num,1,0,AOAspread);

    % AoD of rays with disired seperation
    theta = zeros(ray_num,1);
    theta0 = 0/180*pi;
    theta = theta0 + laprnd(ray_num,1,0,AODspread);

    % Gain
    g_cmplx = exp(1j*rand*2*pi)/ray_num;
    g = g_cmplx;

    % Pre-compute some vectors/matrices in FIM
    for rayindex = 1:ray_num

        % Spatial response and its derivative over phi
        arx(:,rayindex) = exp(1j * pi * (0:Nr-1)' * sin(phi(rayindex)))/sqrt(Nr);
        Darx(:,rayindex) = 1j * pi * (0:Nr-1)' * cos(phi(rayindex));
        drx(:,rayindex) = Darx(:,rayindex).*arx(:,rayindex);

        atx(:,rayindex) = exp(1j * pi * (0:Nt-1)' * sin(theta(rayindex)))/sqrt(Nt);
        Datx(:,rayindex) = 1j * pi * (0:Nt-1)' * cos(theta(rayindex));
        dtx(:,rayindex) = Datx(:,rayindex).*atx(:,rayindex);
    end


    % Evaluation of FIM with multiple rays

    % Pre-compute some equations, vecx, vecalpha, vecbeta
    % vecx is derivative of phi_0 in f(r,eta)
    vaa = zeros(M,1);
    vad = zeros(M,1);
    vda = zeros(M,1);

    for rr=1:ray_num
        vaa = vaa + (W'*arx(:,rr)).* conj(F'*atx(:,rr));
        vad = vad + (W'*arx(:,rr)).* conj(F'*dtx(:,rr));
        vda = vda + (W'*drx(:,rr)).* conj(F'*atx(:,rr));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluation of J_{0,0} part, a 4 by 4 matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J_00 = zeros(4,4);

    % Partial of phi_0 and phi_0              
    J_00(1,1) = real((g*vda)' * (g*vda));

    % Partial of phi_0 and theta_0
    J_00(1,2) = real((g*vda)' * (g*vad));
    J_00(2,1) = J_00(1,2);

    % Partial of phi_0 and alpha
    J_00(1,3) = real((g*vda)' * (1*vaa));
    J_00(3,1) = J_00(1,3);

    % Partial of phi_0 and beta
    J_00(1,4) = real((g*vda)' * (1j*vaa));
    J_00(4,1) = J_00(1,4);

    % Partial theta_0 and theta_0
    J_00(2,2) = real((g*vad)' * (g*vad));

    % Partial theta_0 and alpha
    J_00(2,3) = real((g*vad)' * (1*vaa));
    J_00(3,2) = J_00(2,3);

    % Partial theta_0 and beta
    J_00(2,4) = real((g*vad)' * (1j*vaa));
    J_00(4,2) = J_00(2,4);

    % Partial of alpha and alpha
    J_00(3,3) = real((1*vaa)' * (1*vaa));

    % Partial of alpha and beta
    J_00(3,4) = real((1*vaa)' * (1j*vaa));
    J_00(4,3) = J_00(3,4);

    % Partial of beta and beta
    J_00(4,4) = real((1j*vaa)' * (1j*vaa));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate J_{0,Deltaphi} matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J_0D = zeros(4, ray_num*2);
    for rr=1:ray_num

        vdar = W'*drx(:,rr).* conj(F'*atx(:,rr));
        vadr = W'*arx(:,rr).* conj(F'*dtx(:,rr));

        phir_index = (rr-1)*2+1;
        thetar_index = rr*2;

        % partial phi_0 and phi_r
        J_0D(1,phir_index) = real((g*vda)' * (g*vdar));
        % partial phi_0 and theta_r
        J_0D(1,thetar_index) = real((g*vda)' * (g*vadr));

        % partial theta_0 and phi_r
        J_0D(2,phir_index) = real((g*vad)' * (g*vdar));
        % partial theta_0 and phi_r
        J_0D(2,thetar_index) = real((g*vad)' * (g*vadr));

        % partial alpha and phi_r
        J_0D(3,phir_index) = real((1*vaa)' * (g*vdar));
        % partial alpha and theta_r
        J_0D(3,thetar_index) = real((1*vaa)' * (g*vadr));

        % partial beta and phi_r
        J_0D(4,phir_index) = real((1j*vaa)' * (g*vdar));
        % partial beta and theta_r
        J_0D(4,thetar_index) = real((1j*vaa)' * (g*vadr));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % J(Deltaph,0) is transpose matrix of J_{0,Deltaphi}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J_D0 = J_0D.';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate J{Deltaphi,Deltaphi} matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J_DD = zeros(ray_num*2, ray_num*2);
    for r1 = 1:ray_num
        phir1_index = (r1-1)*2+1;
        thetar1_index = r1*2;
        vadr1 = W'*arx(:,r1).* conj(F'*dtx(:,r1));
        vdar1 = W'*drx(:,r1).* conj(F'*atx(:,r1));
        for r2 = 1:ray_num
            phir2_index = (r2-1)*2+1;
            thetar2_index = r2*2;
            vadr2 = W'*arx(:,r2).* conj(F'*dtx(:,r2));
            vdar2 = W'*drx(:,r2).* conj(F'*atx(:,r2));

            % partial phi_r1 and phi_r2
            J_DD(phir1_index,phir2_index) = real((g*vdar1)'* (g*vdar2));

            % partial phi_r1 and theta_r2
            J_DD(phir1_index,thetar2_index) = real((g*vdar1)'* (g*vadr2));

            % partial theta_r1 and phi_r2
            J_DD(thetar1_index,phir2_index) = real((g*vadr1)'* (g*vdar2));

            % partial theta_r1 and theta_r2
            J_DD(thetar1_index,thetar2_index) = real((g*vadr1)'* (g*vadr2));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate J matrix for apriori distribution of dphi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J_p = zeros(4+ray_num*2);
    % A-priori of angular spread
    J_p(5:end,5:end) = diag(ones(2*ray_num,1)*sqrt(2)/AOAspread);
    
    % For loop for SNR
    for ss = 1:SNR_num
        
        % SNR
        sigman2 = 10^(-SNR_range(ss)/10);
        
        % Evaluate FIM
        J = 2/sigman2*[J_00, J_0D; J_D0, J_DD] + J_p;
        
        % RMSE evaluation from CRLB perspective
        % CRLB of first ray when there are multiple rays
        temp = inv(J);
        CRLB_multiple(ss,MCindex) = sqrt(temp(1,1))*(1/pi*180);
        
        % Evaluation of FIM with single rays
        temp = inv(J(1:4,1:4));
        % RMSE evaluation from CRLB perspective
        CRLB_single(ss,MCindex) = sqrt(temp(1,1))*(1/pi*180);

    end
end
%% Plot CRLB of angle est./tracking
figure
semilogy(SNR_range,mean(CRLB_single,2));hold on
semilogy(SNR_range,mean(CRLB_multiple,2));hold on
grid on
legend('Single Ray','Multple Rays')

% figure
% semilogy(SNR_range,mean(CRLB_rest1,2));hold on
% semilogy(SNR_range,mean(CRLB_rest1,2));hold on
% grid on
% legend('CRLB Ray 1','CRLB Ray 2')
