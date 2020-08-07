function [ A, B, C ] = Get_ABCmathcal( X, time, userindex, h_channel, f_channel, sigma, Method )
%GET_ABCMATHCAL Summary of this function goes here
%   Detailed explanation goes here

% X = X1 or X2

gamma_DL = Get_gammaDL(X, userindex, h_channel, f_channel, sigma, Method);

A = real(2*log(1+gamma_DL)/time + gamma_DL/(time*(gamma_DL+1)));

B = real((gamma_DL.^2)/(time*(gamma_DL+1)));

C = real(log(1+gamma_DL)/(time.^2));

end

function [gamma_DL] = Get_gammaDL( X, userindex, h_channel, f_channel, sigma, Method )

% This function is to calculate gamma_DL as in (24) (29)

[varphi, ~] = Get_varphi(X, userindex, h_channel, f_channel, sigma, Method); % interference plus noise in (6a) (6b)

W = X{1,1};

gamma_DL = (real(h_channel'*W(:,userindex)))^2/varphi;

end

