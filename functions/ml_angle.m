function angle_error = ml_angle(y, true_rayAOA, true_rayAOD, F, W, cluster_num, Nt, Nr)
    low_limit = (true_rayAOA - pi/4)/pi*180;
    up_limit = (true_rayAOA + pi/4)/pi*180;
    resolution = low_limit:3:up_limit;
    resolution'
    
    angle_estimates = zeros(1, length(resolution));
    for cluster_index = 1:cluster_num
        atx = exp(1j*(0:Nt-1)'*pi*sin(true_rayAOD(1,cluster_index)))/sqrt(Nt);
        for angle_index = 1:length(resolution)
          arx = exp(1j*(0:Nr-1)'*pi*sin(resolution(angle_index)/180*pi))/sqrt(Nr);
%             arx = exp(1j*(0:Nr-1)'*pi*sin(true_rayAOA))/sqrt(Nr);

            x = W'*arx*atx'*F*ones(4,1);
%             angle_estimates = (abs(y'*x))^2/norm(x)^2;
            angle_estimates(1, angle_index) = (abs(y'*x))^2/norm(x)^2;
        end
        [~,index] = max(angle_estimates);
        angle_est = resolution(index); %estimated w=(w1,w2)
        true_rayAOA = true_rayAOA/pi*180;
        angle_error = abs(true_rayAOA - angle_est);
    end
end