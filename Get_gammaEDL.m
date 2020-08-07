function [ gamma_EDL, psi_tilde_m ] = Get_gammaEDL( X, time, userindex, Hm, sigma, Lambda, zone )
%GET_GAMMAEDL Summary of this function goes here
%   Detailed explanation goes here

N_em = size(Hm,2);
w = X{1,zone}{1,1}(:,userindex);

time_vec = mod(eye(2)+(zone-1),2)*[1; (time-1)];
% time_vec = [1; 1/(time-1)];
% Fact = [1 (time-1)];

psi_tilde_m = Lambda*time_vec + time*N_em*sigma^2 - (norm(Hm'*w))^2;
% psi_tilde_m = Fact(zone)*(Lambda*time_vec + time/(time-1)*N_em*sigma^2) - (norm(Hm'*w))^2;
 
gamma_EDL = (norm(Hm'*w))^2/psi_tilde_m;

end

