function [ A_tilde, B_tilde, C_tilde, Omega ] = Get_ABCtildemathcal( X, time, userindex, G_channel, G_SI, sigma_SI, sigma_noise )
%GET_ABCTILDEMATHCAL Summary of this function goes here
%   Detailed explanation goes here

% This function is to get parameters for (33) (34)
% X is X1 or X2, time is alpha or beta

[gamma_UL, Omega] = Get_gammaUL(X, userindex, G_channel, G_SI, sigma_SI, sigma_noise);

A_tilde = real((2*log(1+gamma_UL) - gamma_UL)/time);

B_tilde = real((2*gamma_UL)/(X{1,3}(userindex)*time));

C_tilde = real(log(1+gamma_UL)/(time.^2));

end


function [ gamma_UL, Omega ] = Get_gammaUL( X, userindex, G_channel, G_SI, sigma_SI, sigma_noise )

Phi = Get_Phi(X, userindex, G_channel(:,2:end), G_SI, sigma_SI, sigma_noise );

gamma_UL = (X{1,3}(userindex))^2*(G_channel(:,1))'*inv(Phi)*G_channel(:,1);

Omega = inv(Phi) - inv((X{1,3}(userindex))^2*G_channel(:,1)*(G_channel(:,1))' + Phi);

end

function [ Phi ] = Get_Phi(X, userindex, G_channel, G_SI, sigma_SI, sigma_noise )

temp = 0;


for iULU = 1:1:size(G_channel,2)
    
    temp = temp + (X{1,3}(userindex + iULU))^2*G_channel(:,iULU)*(G_channel(:,iULU))';
end

temp = temp + sigma_SI*G_SI'*(X{1,1}*(X{1,1})'+ X{1,2}*(X{1,2})')*G_SI;

Phi = temp + sigma_noise^2*eye(size(temp,1));


end