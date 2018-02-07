%-------------------------------------
% Script control parameter
%-------------------------------------
clear;clc;
rng(2); %random seed

%-------------------------------------
% System Parameters
%-------------------------------------
ray_num = 2;
Nr = 16;
M = 5;

for MCindex = 1:100
    
    % Receiver beamformer: 1) quasi-omni beam from random steering mtx; 2)
    % directional beam from angle steering vector
    W = (randi(2,Nr,M)*2-3) + 1j * (randi(2,Nr,M)*2-3);
    % w = exp(1j * pi * (0:Nr-1)' * sin(0));
    % W = repmat(w,1,M);
    
    % Evaluate CRLB when the second ray is phi apart
    dphi_num = 100;
    dphi_range = linspace(2,15,dphi_num);

    for pp = 1:dphi_num
        phi = zeros(ray_num,1);
        phi(1) = 0/180*pi;
        phi(2) = phi(1) + dphi_range(pp)/180*pi;
        
        % Pre-compute some vectors/matrices in FIM
        for rayindex = 1:ray_num
            
            % Complex gain
            g_cmplx(rayindex) = (randn+1j*randn)/sqrt(ray_num);

            % Spatial response and its derivative over phi
            arx(:,rayindex) = exp(1j * pi * (0:Nr-1)' * sin(phi(rayindex)));
            Darx(:,rayindex) = 1j * pi * (0:Nr-1)' * cos(phi(rayindex));
        end
        
        % Evaluation of FIM
        for r1 = 1:ray_num
            for r2 = 1:ray_num
                % Partial of phi_r1 and phi_r2
                J_sub(1,1) = real((g_cmplx(r1)*W'*diag(Darx(:,r1))*arx(:,r1))'...
                                *(g_cmplx(r2)*W'*diag(Darx(:,r2))*arx(:,r2)));

                % Partial of phi_r1 and alpha_r2
                J_sub(1,2) = real((g_cmplx(r1)*W'*diag(Darx(:,r1))*arx(:,r1))'...
                                *(W'*arx(:,r2)));

                % Partial of phi_r1 and beta_r2
                J_sub(1,3) = real((g_cmplx(r1)*W'*diag(Darx(:,r1))*arx(:,r1))'...
                                *(1j*W'*arx(:,r2)));

                % Partial of alpha_r1 and phi_r2
                J_sub(2,1) = real((W'*arx(:,r1))'...
                                *(g_cmplx(r2)*W'*diag(Darx(:,r2))*arx(:,r2)));

                % Partial of alpha_r1 and alpha_r2
                J_sub(2,2) = real((W'*arx(:,r1))'...
                                *(W'*arx(:,r2)));

                % Partial of alpha_r1 and beta_r2
                J_sub(2,3) = real((W'*arx(:,r1))'...
                                *(1j*W'*arx(:,r2)));

                % Partial of beta_r1 and phi_r2
                J_sub(3,1) = real((1j*W'*arx(:,r1))'...
                                *(g_cmplx(r2)*W'*diag(Darx(:,r2))*arx(:,r2)));

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

        % CRLB of first ray when there is only one ray
        xindex = 1:3;
        yindex = 1:3;
        temp = inv(J(xindex,yindex));
        % RMSE evaluation from CRLB perspective
        CRLB_single(pp,MCindex) = sqrt(temp(1,1))*(1/pi*180);
%         CRLB_single(pp,MCindex) = sqrt(1/(J(1,1) -...
%         J(1,2:3)*inv(J(2:3,2:3))*J(2:3,1)))*(1/pi*180)

        % RMSE evaluation from CRLB perspective
        % CRLB of first ray when there are multiple rays
        temp = inv(J);
        CRLB_multiple(pp,MCindex) = sqrt(temp(1,1))*(1/pi*180);
    end
end
%%
figure
semilogy(dphi_range,mean(CRLB_single,2));hold on
semilogy(dphi_range,mean(CRLB_multiple,2));hold on
grid on

