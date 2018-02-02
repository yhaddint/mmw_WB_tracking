function [ raygain, raydelay, ray_AOA_azim, ray_AOD_azim ] =...
    get_chan_parameter_nogeo( print_stat,...
                              cluster_num,...
                              ray_num,...
                              sigma_delay_spread,...
                              centroid_AOA,...
                              sigma_AOA_spread,...
                              centroid_AOD,...
                              sigma_AOD_spread )
%GET_CHAN_PARAMETER_NOGEO Summary of this function goes here
%   Generate channel parameter using given statistics
%   [ raygain, raydelay, ray_AOA_azim, ray_AOD_azim ] = ...
%           get_chan_parameter_nogeo( print_stat, cluster_num, ray_num, sigma_delay_spread,...
%                                      centroid_AOA, sigma_AOA_spread,...
%                                      centroid_AOD, sigma_AOD_spread )
%   IP: print_stat, indicator for printing channel statistics
%   IP: cluster_num, number of MPC cluster
%   IP: ray_num, number of rays in each clusters
%   IP: sigma_delay_spread, desired intra-cluster delay spread (std with
%       unit second)
%   IP: centroid_AOA, desired centroid of AOA of clusters (with unit rad), can
%       be 'random' if not to be fixed
%   IP: sigma_AOA_spread, desired intra-cluster AOA spread (std var with unit rad)
%   IP: centroid_AOD, desired centroid of AOD of clusters (with unit rad), can
%       be 'random' if not to be fixed
%   IP: sigma_AOD_spread, desired intra-cluster AOD spread (std var with unit rad)
%   OP: raygain, cluster_num by ray_num matrix with raygain for each ray
%   OP: raydelay, cluster_num by ray_num matrix with delay for each ray
%   OP: ray_AOA_azim, cluster_num by ray_num matrix with AOA for each ray
%   OP: ray_AOD_azim, cluster_num by ray_num matrix with AOD for each ray


if cluster_num>2
    fprintf('This function only support up to 2 multipath clusters\n');
end

if max(sigma_delay_spread)>1e-6
    fprintf('Delay spread has beter not to exceed 1e-6 second\n');
end

cluster_delay = [300e-9,250e-9]; % Mean delay of two multipath clusters

% Zero initializations
raygain = zeros(cluster_num, ray_num);
raydelay = zeros(cluster_num, ray_num);
rayAOA = zeros(cluster_num, ray_num);
rayAOD = zeros(ray_num, cluster_num);


for cluster_index = 1:cluster_num

    % Randomly generate chan. parameters (AOD in Azimuth)
    if centroid_AOD == 'random'
        ray_AOD_azim_mean = rand * 2 * pi / 3 - pi/3;
    else
        ray_AOD_azim_mean = centroid_AOD(cluster_index);
    end
    angle_spread = laprnd(1, ray_num, 0, sigma_AOD_spread);
    ray_AOD_azim(cluster_index,:) = ray_AOD_azim_mean + angle_spread;

%     % Randomly generate chan. parameters (AOD in Elevation)
%     ray_AOD_elev_mean = rand * pi / 3 - pi/6;
%     angle_spread = laprnd(1, ray_num, 0, sigma_AOD_spread);
%     ray_AOD_elev(cluster_index,:) = ray_AOD_elev_mean + angle_spread;

    % Randomly generate chan. parameters (AOA in Azimuth)
    if centroid_AOA == 'random'
        ray_AOA_azim_mean = rand * 2 * pi / 3 - pi/3;
    else
        ray_AOA_azim_mean = centroid_AOA(cluster_index);
    end
    angle_spread = laprnd(1, ray_num, 0, sigma_AOA_spread);
    ray_AOA_azim(cluster_index,:) = ray_AOA_azim_mean + angle_spread;

%     % Randomly generate chan. parameters (AOA in Elevation)
%     ray_AOA_elev_mean = rand * pi / 3 - pi/6;
%     angle_spread = laprnd(1, ray_num, 0, sigma_AOA_spread);
%     ray_AOA_elev(cluster_index,:) = ray_AOA_elev_mean + angle_spread;


    % Unit gain for each ray
    raygain(cluster_index,:) = ones(1,ray_num);
    
    % Delay of each ray, Gaussian distributed
    raydelay(cluster_index,:) = cluster_delay(cluster_index) + randn(1,ray_num)*sigma_delay_spread;
    
    if print_stat
        fprintf('Cluster %d:\n',cluster_index);
        fprintf('DS Mean:    %4.2f ns\n', mean(raydelay(cluster_index,:)*1e9));
        fprintf('DS Std Dev: %4.2f ns\n', sqrt(var(raydelay(cluster_index,:)*1e9)));
        fprintf('AOAS Mean:    %4.2f degree\n', mean(ray_AOA_azim(cluster_index,:)/pi*180));
        fprintf('AOAS Std Dev: %4.2f degree\n', sqrt(var(ray_AOA_azim(cluster_index,:)/pi*180)));
        fprintf('AODS Mean:    %4.2f degree\n', mean(ray_AOD_azim(cluster_index,:)/pi*180));
        fprintf('AODS Std Dev: %4.2f degree\n', sqrt(var(ray_AOD_azim(cluster_index,:)/pi*180)));
    end
end

end

