function [ varphi, varphi_vec ] = Get_varphi( X, userindex, h_channel, f_channel, sigma, Method )
%GET_VARPHI Summary of this function goes here
%   Detailed explanation goes here

% varphi = norm(varphi_vec)^2 - Interference plus noise in (6)


temp = 0;
temp_vec = [];

W = X{1,1};
V = X{1,2};
rho = X{1,3};
K2 = floor(size(W,2)/2);

for iDLU = 1:1:size(W,2)
    
    if ((iDLU==userindex) || ((Method==2) && (userindex<=K2) && (iDLU==(K2+userindex)) ) ) % Method=2 for NOMA
        continue
    end
    
    temp = temp + (abs(h_channel'*W(:,iDLU)))^2;
    temp_vec = [temp_vec, (h_channel'*W(:,iDLU))];
    
end

temp = temp + (norm(h_channel'*V))^2;
temp_vec = [temp_vec, (h_channel'*V)];
if (Method~=3)
    for iULU = 1:1:length(f_channel)
        temp = temp + (rho(iULU))^2*(abs(f_channel(iULU)))^2;
        temp_vec = [temp_vec, (rho(iULU)*f_channel(iULU))];
    end
end

varphi = temp + sigma^2;
varphi_vec = [temp_vec, sigma];

end

