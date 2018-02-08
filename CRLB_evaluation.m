%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;
rng(2); %random seed
% load probe_BF
%-------------------------------------
% System Parameters
%-------------------------------------
ray_num = 2; % Num of rays in a cluster
Nr = 32; % Number of antenna in Rx
M = 4; % Length of training
MCtimes = 50; % Num of Monte Carlo Sim.


% For loop for Monte Carlo Simulations
for MCindex = 1:MCtimes
    
    % Receiver beamformer: 1) quasi-omni beam from random steering mtx; 2)
    % directional beam from angle steering vector
    probe_BF = (randi(2,Nr,M)*2-3) + 1j * (randi(2,Nr,M)*2-3);
    W = probe_BF./norm(probe_BF,'fro');
%     W = probe_BF(:,1:4)./norm(probe_BF(:,1:4),'fro');
    
    % Evaluate CRLB when the second ray is phi apart
    dphi_num = 100;
    dphi_range = linspace(2,15,dphi_num);
    
    % For loop for each two-ray seperation
    for pp = 1:dphi_num
        
        % AoA of rays with disired seperation
        phi = zeros(ray_num,1);
        phi(1) = 30/180*pi;
%         randompm = randi(2,ray_num-1,1)*2-3;
        phi(2) = phi(1) + dphi_range(pp)/180*pi;
        phi(3) = phi(1) + 2*dphi_range(pp)/180*pi;
        
        % Pre-compute some vectors/matrices in FIM
        for rayindex = 1:ray_num
            
            % Complex gain
%             g_cmplx(rayindex) = (randn + 1j*randn)/sqrt(ray_num);
            g_cmplx(rayindex) = exp(1j*rand*2*pi);

            % Spatial response and its derivative over phi
            arx(:,rayindex) = exp(1j * pi * (0:Nr-1)' * sin(phi(rayindex)));
            Darx(:,rayindex) = 1j * pi * (0:Nr-1)' * cos(phi(rayindex));
        end
        
        % Evaluation of FIM with multiple rays
        vecx = zeros(M,1);
        for rr=1:ray_num
            vecx = vecx + g_cmplx(rr)*W'*(Darx(:,rr).*arx(:,rr));
        end
        for r1 = 1:ray_num
            for r2 = 1:ray_num
                
                % Partial of phi_r1 and phi_r2
                if r1 == 1
                    vec1 = vecx;
                else
                    vec1 = g_cmplx(r1)*W'*diag(Darx(:,r1))*arx(:,r1);
                end
                if r2 == 1
                    vec2 = vecx;
                else
                    vec2 = g_cmplx(r2)*W'*diag(Darx(:,r2))*arx(:,r2);
                end
                J_sub(1,1) = real(vec1' * vec2);

                % Partial of phi_r1 and alpha_r2
                if r1 == 1
                    vec1 = vecx;
                else
                    vec1 = g_cmplx(r1)*W'*diag(Darx(:,r1))*arx(:,r1);
                end
                vec2 = (W'*arx(:,r2));
                J_sub(1,2) = real(vec1' * vec2);

                % Partial of phi_r1 and beta_r2
                if r1 == 1
                    vec1 = vecx;
                else
                    vec1 = g_cmplx(r1)*W'*diag(Darx(:,r1))*arx(:,r1);
                end
                vec2 = 1j*W'*arx(:,r2);
                J_sub(1,3) = real(vec1'*vec2);

                % Partial of alpha_r1 and phi_r2
                vec1 = (W'*arx(:,r1));
                if r2 == 1
                    vec2 = vecx;
                else
                    vec2 = g_cmplx(r2)*W'*diag(Darx(:,r2))*arx(:,r2);
                end
                J_sub(2,1) = real(vec1'*vec2);

                % Partial of alpha_r1 and alpha_r2
                J_sub(2,2) = real((W'*arx(:,r1))'...
                                *(W'*arx(:,r2)));

                % Partial of alpha_r1 and beta_r2
                J_sub(2,3) = real((W'*arx(:,r1))'...
                                *(1j*W'*arx(:,r2)));

                % Partial of beta_r1 and phi_r2
                vec1 = (1j*W'*arx(:,r1));
                if r2 == 1
                    vec2 = vecx;
                else
                    vec2 = g_cmplx(r2)*W'*diag(Darx(:,r2))*arx(:,r2);
                end
                J_sub(3,1) = real(vec1'* vec2);

                % Partial of beta_r1 and alpha_r2
                J_sub(3,2) = real((1j*W'*arx(:,r1))'...
                                *(W'*arx(:,r2)));

                % Partial of beta_r1 and beta_r2
                J_sub(3,3) = real((1j*W'*arx(:,r1))'...
                                *(1j*W'*arx(:,r2)));

                % indecies that J_sub should go into J matrix
                xindex = (r1-1)*3+1:r1*3;
                yindex = (r2-1)*3+1:r2*3;

                % Fill into J matrix
                J(xindex,yindex) = J_sub;
            end
        end
        
        % RMSE evaluation from CRLB perspective
        % CRLB of first ray when there are multiple rays
        temp = inv(J);
        CRLB_multiple(pp,MCindex) = sqrt(temp(1,1))*(1/pi*180);
        
        % Evaluation of FIM with single rays
        for r1 = 1:1
            for r2 = 1:1
                
                % Partial of phi_r1 and phi_r2
                vec1 = g_cmplx(1)*W'*diag(Darx(:,1))*arx(:,1);
                vec2 = g_cmplx(1)*W'*diag(Darx(:,1))*arx(:,1);
                J_sub(1,1) = real(vec1' * vec2);

                % Partial of phi_r1 and alpha_r2
                vec1 = g_cmplx(1)*W'*diag(Darx(:,1))*arx(:,1);
                vec2 = (W'*arx(:,1));
                J_sub(1,2) = real(vec1' * vec2);

                % Partial of phi_r1 and beta_r2
                vec1 = g_cmplx(1)*W'*diag(Darx(:,1))*arx(:,1);
                vec2 = 1j*W'*arx(:,1);
                J_sub(1,3) = real(vec1'*vec2);

                % Partial of alpha_r1 and phi_r2
                vec1 = (W'*arx(:,1));
                vec2 = g_cmplx(1)*W'*diag(Darx(:,1))*arx(:,1);
                J_sub(2,1) = real(vec1'*vec2);

                % Partial of alpha_r1 and alpha_r2
                J_sub(2,2) = real((W'*arx(:,1))'...
                                *(W'*arx(:,1)));

                % Partial of alpha_r1 and beta_r2
                J_sub(2,3) = real((W'*arx(:,1))'...
                                *(1j*W'*arx(:,1)));

                % Partial of beta_r1 and phi_r2
                vec1 = 1j*W'*arx(:,1);
                vec2 = g_cmplx(1)*W'*diag(Darx(:,1))*arx(:,1);
                J_sub(3,1) = real(vec1'* vec2);

                % Partial of beta_r1 and alpha_r2
                J_sub(3,2) = real((1j*W'*arx(:,1))'...
                                *(W'*arx(:,1)));

                % Partial of beta_r1 and beta_r2
                J_sub(3,3) = real((1j*W'*arx(:,1))'...
                                *(1j*W'*arx(:,1)));

                % Fill into J matrix
                J = J_sub;
            end
        end

        % CRLB of first ray when there is only one ray
        xindex = 1:3;
        yindex = 1:3;
        temp = inv(J(xindex,yindex));
        % RMSE evaluation from CRLB perspective
        CRLB_single(pp,MCindex) = sqrt(temp(1,1))*(1/pi*180);
%         CRLB_single(pp,MCindex) = sqrt(1/(J(1,1) - J(1,2:3)*inv(J(2:3,2:3))*J(2:3,1)))*(1/pi*180);


    end
end
%% Plot CRLB of angle est./tracking
figure
semilogy(dphi_range,mean(CRLB_single,2));hold on
semilogy(dphi_range,mean(CRLB_multiple,2));hold on
grid on

