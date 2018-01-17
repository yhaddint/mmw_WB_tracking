function [tau_est, deltaAOA_est,  AOA_est] = LMS_tracking_v2( alpha, tau, phi, phi_steer, theta, theta_steer, h_now , Nt, Nr, Nfft, ray_num, stepsize)
%LMS_TRACKING Summary of this function goes here
%   Detailed explanation goes here
 % Alternative estimation parameters
        
        alpha_est = alpha;
        %---------------------------------------------
        % tau estimation with previous alpha and phi
        %---------------------------------------------
        for kk=1:Nfft
            bigPhi = sin(phi)-sin(phi_steer);
            bigTheta = -sin(theta)+sin(theta_steer);
            BigTau(kk,:) = (-1j*2*pi/(Nfft*1e-9)*kk)...
            .*exp(-1j*2*pi*kk*tau/(1e-9*Nfft))...
            .*alpha...    
            .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
            .*Nt...
            ./(Nt*Nr);
            
            BigTau_const(kk,:) = exp(-1j*2*pi*kk*tau/(1e-9*Nfft))...
            .*alpha...    
            .*(1-exp(1j*pi*Nr*bigPhi))./(1-exp(1j*pi*bigPhi))...
            .*Nt...
            ./(Nt*Nr);
        end
        deltatau_est = pinv([real(BigTau);imag(BigTau)])...
            *[real(h_now-sum(BigTau_const,2));imag(h_now-sum(BigTau_const,2))];
        tau_est = tau + deltatau_est.';
        
        
        %---------------------------------------------
        % phi estimation using two assisted receiving beams
        %---------------------------------------------
        for kk=1:Nfft
            bigPhi_bf1 = sin(phi)-sin(phi_steer);
            
            bigTheta = -sin(theta)+sin(theta_steer);
            BigPhi_const_bf1(kk,:) = alpha(1,:)...
                .*exp(-1j*2*pi*kk*tau_est(1,:)/(1e-9*Nfft))...
                .*(1-exp(1j*pi*Nr*bigPhi_bf1))./(1-exp(1j*pi*bigPhi_bf1))...
                .*Nt...
                ./(Nt*Nr);
            
            
            for ray_index=1:ray_num
                BigPHI_bf1(kk,ray_index) = alpha(ray_index)...
                .*exp(-1j*2*pi*kk*tau_est(ray_index)/(1e-9*Nfft))...
                .*(((1-exp(1j*pi*bigPhi_bf1(ray_index)))*(-exp(1j*pi*Nr*bigPhi_bf1(ray_index)))*(1j*pi*Nr*cos(phi(ray_index))))...
                  -((1-exp(1j*pi*Nr*bigPhi_bf1(ray_index)))*(-exp(1j*pi*bigPhi_bf1(ray_index)))*(1j*pi*cos(phi(ray_index)))))...
                ./((1-exp(-1j*pi*bigPhi_bf1(ray_index))).^2)...
                .*Nt...
                ./(Nt*Nr);
            
            end
        end
        
        y1_LS_comp = h_now - sum(BigPhi_const_bf1,2);
        A1_LS_comp  = BigPHI_bf1;
        
        y_LS = [real(y1_LS_comp);imag(y1_LS_comp)];
        A_LS = [real(A1_LS_comp);imag(A1_LS_comp)];
        
        deltaAOA_est = sum(A_LS,2)'*y_LS;
        AOA_est = phi + stepsize * deltaAOA_est.';

        

end

