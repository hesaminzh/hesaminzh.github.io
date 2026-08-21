function [user_loc, user_azimuths, d_BR, d_RU, d_BU] = generate_user_location(user_range, user_center, RIS_loc, BS_loc, K)
%Generate the user locations randomly 
%  [user_loc, user_angle, d_RU, d_BU] = generate_user_location(para)
%Inputs:
%   para: structure of the initial parameters
%Outputs:
%   user_loc: locations of users in km
%   user_agnle: angle of users from RIS
%   d_RU: distance between RIS and users
%   d_BU: distance between BS and RIS
%Date: 30/05/2021
%Author: Zhaolin Wang
% modified by Hessam Nezhadmohammad



d_CU = rand(K,1) * (user_range(2) - user_range(1))...
    + user_range(1); % distance from user center to user
user_angle_center = rand(K,1) * (2 *pi); % angle of directions from RIS
user_loc = d_CU.*exp(1i*user_angle_center);
user_loc = [real(user_loc), imag(user_loc)] + user_center;
% user_loc = [206.5590, 35.2627; 201.9195, 35.1499; 205.8163, 28.9895; 207.0316, 27.9245];

relative_user_loc = user_loc - RIS_loc;
d_BR = sqrt(sum((RIS_loc - BS_loc).^2));
d_RU = sqrt(relative_user_loc(:,1).^2 + relative_user_loc(:,2).^2);
d_BU = sqrt(user_loc(:,1).^2 + user_loc(:,2).^2);

%user_angle = atan(relative_user_loc(:,2) ./ relative_user_loc(:,1)) - pi/4;
%user_angle = atan2(relative_user_loc(:,2), relative_user_loc(:,1));
user_azimuths = atan2( user_loc(:,2) - RIS_loc(2), user_loc(:,1) - RIS_loc(1) );

end

