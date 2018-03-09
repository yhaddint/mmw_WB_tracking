function [raydelay, rayAOA, rayAOD ] = get_multipath_v2(loc_bs, loc_ue, loc_cluster_total,...
                                         cluster_num, ray_num, print_stat )
%GET_MULTIPATH Summary of this function goes here
%   Detailed explanation goes here

    % adjusted the orientation of Tx array

    speed_c = 3e8;
    for cluster_index = 1:cluster_num
       for ray_index = 1:ray_num
           loc_cluster = loc_cluster_total((cluster_index-1)*ray_num+1:cluster_index*ray_num,:);
            raydelay(cluster_index,ray_index) = (norm(loc_cluster(ray_index,:)-loc_bs)...
                        +norm(loc_cluster(ray_index,:)-loc_ue))/speed_c;
           temp = loc_cluster(ray_index,:) - loc_ue;
            rayAOA(cluster_index,ray_index) = angle(temp(1)+1j*temp(2))-30/180*pi;
            temp = loc_cluster(ray_index,:) - loc_bs;
%             rayAOD(cluster_index,ray_index) = angle(-temp(1)-1j*temp(2));
            rayAOD(cluster_index,ray_index) = 0; % unrealistic case with 0
       end
        if print_stat
            clc
            fprintf('Cluster %d:\n',cluster_index);
            fprintf('DS Mean:    %4.2f ns\n', mean(raydelay(cluster_index,:)*1e9));
            fprintf('DS Std Dev: %4.2f ns\n', sqrt(var(raydelay(cluster_index,:)*1e9)));
            fprintf('AOAS Mean:    %4.2f degree\n', mean(rayAOA(cluster_index,:)/pi*180));
            fprintf('AOAS Std Dev: %4.2f degree\n', sqrt(var(rayAOA(cluster_index,:)/pi*180)));
            fprintf('AODS Mean:    %4.2f degree\n', mean(rayAOD(cluster_index,:)/pi*180));
            fprintf('AODS Std Dev: %4.2f degree\n', sqrt(var(rayAOD(cluster_index,:)/pi*180)));
        end
    end

end

