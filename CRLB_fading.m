%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;
rng(3); %random seed
% load probe_BF
%-------------------------------------
% System Parameters
%-------------------------------------
ray_num = 15; % Num of rays in a cluster
Nr = 8; % Number of antenna in Rx
Nt = 32;
M = 64; % Length of training
MCtimes = 20; % Num of Monte Carlo Sim.
AOAspread2 = (5/180*pi)^2;
AOAspread = sqrt(AOAspread2);
AODspread2 = (5/180*pi)^2;
AODspread = sqrt(AODspread2);
SNR_num = 50;
SNR_range = linspace(-15,30,SNR_num);

% For loop for Monte Carlo Simulations (realization of g, angle spread, and beamformer)
for MCindex = 1:MCtimes
    
    clc
    fprintf('iteration %d:\n',MCindex);
    
    % Receiver beamformer: 1) quasi-omni beam from random steering mtx; 2)
    % directional beam from angle steering vector
    probe_Rx_BF = (randi(2,Nr,M)*2-3) + 1j * (randi(2,Nr,M)*2-3);
    W = probe_Rx_BF./norm(probe_Rx_BF,'fro')*sqrt(M);
    
    probe_Tx_BF = (randi(2,Nt,M)*2-3) + 1j * (randi(2,Nt,M)*2-3);
    F = probe_Tx_BF./norm(probe_Tx_BF,'fro')*sqrt(Nt*M);
    
%     probe_Tx_BF = ones(Nt,M);
%     F = probe_Tx_BF./norm(probe_Tx_BF,'fro')*sqrt(Nt*M);   

    % AoA of rays with disired seperation
    phi = zeros(ray_num,1);
    phi0 = 0/180*pi;
    phi = phi0 + randn(ray_num,1) * AOAspread;

    % AoD of rays with disired seperation
    theta = zeros(ray_num,1);
    theta0 = 0/180*pi;
    theta = theta0 + randn(ray_num,1) * AODspread;

%     % Gain
%     g_cmplx = exp(1j*rand(ray_num,1)*2*pi)/sqrt(ray_num);
%     g = g_cmplx;
    % Rotate of ray
    psi = rand(ray_num,1)*2*pi;
    g_ray = exp(1j*psi)/sqrt(ray_num);



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
        vaa = vaa + g_ray(rr) * (W'*arx(:,rr)).* conj(F'*atx(:,rr));
        vad = vad + g_ray(rr) * (W'*arx(:,rr)).* conj(F'*dtx(:,rr));
        vda = vda + g_ray(rr) * (W'*drx(:,rr)).* conj(F'*atx(:,rr));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluation of J_{0,0} part, a 3 by 3 matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J_00 = zeros(3,3);

    % Partial of phi_0 and phi_0              
    J_00(1,1) = real((vda)' * (vda));

    % Partial of phi_0 and theta_0
    J_00(1,2) = real((vda)' * (vad));
    J_00(2,1) = J_00(1,2);

    % Partial of phi_0 and alpha
    J_00(1,3) = real((vda)' * (1*vaa));
    J_00(3,1) = J_00(1,3);

    % Partial theta_0 and theta_0
    J_00(2,2) = real((vad)' * (vad));

    % Partial theta_0 and alpha
    J_00(2,3) = real((vad)' * (1*vaa));
    J_00(3,2) = J_00(2,3);

    % Partial of alpha and alpha
    J_00(3,3) = real((1*vaa)' * (1*vaa));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate J_{0,Deltaphi} matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J_0D = zeros(3, ray_num*3);
    for rr=1:ray_num

        vdar = g_ray(rr) * W'*drx(:,rr).* conj(F'*atx(:,rr));
        vadr = g_ray(rr) * W'*arx(:,rr).* conj(F'*dtx(:,rr));
        vpsir = (1j) * psi(rr) * W'*arx(:,rr).* conj(F'*atx(:,rr));

        phir_index = (rr-1)*3+1;
        thetar_index = (rr-1)*3+2;
        psir_index = rr*3;

        % partial phi_0 and phi_r
        J_0D(1,phir_index) = real((vda)' * (vdar));
        % partial phi_0 and theta_r
        J_0D(1,thetar_index) = real((vda)' * (vadr));
        % partial phi_0 and psi_r
        J_0D(1,psir_index) = real((vda)' * (vpsir));

        % partial theta_0 and phi_r
        J_0D(2,phir_index) = real((vad)' * (vdar));
        % partial theta_0 and phi_r
        J_0D(2,thetar_index) = real((vad)' * (vadr));
        % partial theta_0 and psi_r
        J_0D(2,psir_index) = real((vad)' * (vpsir));

        % partial alpha and phi_r
        J_0D(3,phir_index) = real((1*vaa)' * (vdar));
        % partial alpha and theta_r
        J_0D(3,thetar_index) = real((1*vaa)' * (vadr));
        % partial theta_0 and psi_r
        J_0D(3,psir_index) = real((1*vaa)' * (vpsir));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % J(Deltaph,0) is transpose matrix of J_{0,Deltaphi}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J_D0 = J_0D.';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate J{Deltaphi,Deltaphi} matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J_DD = zeros(ray_num*3, ray_num*3);
    for r1 = 1:ray_num
        phir1_index = (r1-1)*3+1;
        thetar1_index = (r1-1)*3+2;
        psir1_index = r1*3;
        vadr1 = g_ray(r1) * W'*arx(:,r1).* conj(F'*dtx(:,r1));
        vdar1 = g_ray(r1) * W'*drx(:,r1).* conj(F'*atx(:,r1));
        vpsir1 = (1j) * psi(r1) * W'*arx(:,r1).* conj(F'*atx(:,r1));
        for r2 = 1:ray_num
            phir2_index = (r2-1)*3+1;
            thetar2_index = (r2-1)*3+2;
            psir2_index = r2*3;
            vadr2 = g_ray(r2) * W'*arx(:,r2).* conj(F'*dtx(:,r2));
            vdar2 = g_ray(r2) * W'*drx(:,r2).* conj(F'*atx(:,r2));
            vpsir2 = (1j) * psi(r2) * W'*arx(:,r2).* conj(F'*atx(:,r2));

            % partial phi_r1 and phi_r2
            J_DD(phir1_index,phir2_index) = real((vdar1)'* (vdar2));

            % partial phi_r1 and theta_r2
            J_DD(phir1_index,thetar2_index) = real((vdar1)'* (vadr2));
            
            % partial phi_r1 and psi_r2
            J_DD(phir1_index,psir2_index) = real((vdar1)'* (vpsir2));

            % partial theta_r1 and phi_r2
            J_DD(thetar1_index,phir2_index) = real((vadr1)'* (vdar2));

            % partial theta_r1 and theta_r2
            J_DD(thetar1_index,thetar2_index) = real((vadr1)'* (vadr2));
            
            % partial theta_r1 and psi_r2
            J_DD(thetar1_index,psir2_index) = real((vadr1)'* (vpsir2));
            
            % partial psi_r1 and phi_r2
            J_DD(psir1_index,phir2_index) = real((vpsir1)'* (vdar2));

            % partial psi_r1 and theta_r2
            J_DD(psir1_index,thetar2_index) = real((vpsir1)'* (vadr2));
            
            % partial psi_r1 and psi_r2
            J_DD(psir1_index,psir2_index) = real((vpsir1)'* (vpsir2));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate J matrix for apriori distribution of dphi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J_p = zeros(3+ray_num*3);
    % A-priori of angular spread
    for r1 = 4:3:(ray_num*3+3)
        J_p(r1,r1) = 2/AOAspread;
    end
    for r2 = 5:3:(ray_num*3+3)
        if AODspread==0
            J_p(r2,r2) = 0;
        else
            J_p(r2,r2) = 2/AODspread;
        end
    end
    
    % For loop for SNR
    for ss = 1:SNR_num
        
        % SNR
        sigman2 = 10^(-SNR_range(ss)/10);
        
        % Evaluate FIM
        J = 2/sigman2*[J_00, J_0D; J_D0, J_DD] + J_p;
        
        % RMSE evaluation from CRLB perspective
        % CRLB of first ray when there are multiple rays
        dphiindex = 4:3:(ray_num*3+3);
        dpsiindex = 6:3:(ray_num*3+3);
        usefulindex = reshape([dphiindex;dpsiindex],1,ray_num*2);
        temp = inv(J);
%         temp = inv(J([1,3,usefulindex],[1,3,usefulindex]));
        CRLB_multiple(ss,MCindex) = sqrt(temp(1,1))*(1/pi*180);
        
        % Evaluation of FIM with single rays
%         temp = inv(J(1:3,1:3));
        temp = inv([J(1,1),J(1,3);J(3,1),J(3,3)]);
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
xlabel('Point-to-Point SNR [dB]')
ylabel('RMSE of AoA [deg]')

% figure
% semilogy(SNR_range,mean(CRLB_rest1,2));hold on
% semilogy(SNR_range,mean(CRLB_rest1,2));hold on
% grid on
% legend('CRLB Ray 1','CRLB Ray 2')
