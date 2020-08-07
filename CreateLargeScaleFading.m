function [ D_H1k, D_H2l, D_G2p, D_G1q, D_f2_pq_kl, D_f1ql, D_Hm, D_Gm2p, D_Gm1q, positionDownlinkUsers, positionUplinkUsers, positionEavs ] = CreateLargeScaleFading( K, L, P, Q, M, M1, Parameters )
%CREATED Summary of this function goes here
%   Detailed explanation goes here

RadiusOfCell = Parameters(1);
InnerZoneRadius = Parameters(2);
RadiusOfNearestUser = Parameters(3);
StandardDeviation = 10^(Parameters(4)/10);
ploss = Parameters(5);

% L = K + M;


%% large-fading for K downlink users - inner zone

%     Zvector = StandardDeviation*randn(1,K);
    rvector = RadiusOfNearestUser*ones(1,K) + (InnerZoneRadius-RadiusOfNearestUser)*rand(1,K);
    anglevector = 2*pi*rand(1,K);
    
    positionDownlinkUsers = [(rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];
        
    distanceBS_To_DownlinkUsers_Zone1 = sqrt(sum(positionDownlinkUsers'.^2));
        
%     betavector_downlink_Zone1 = (10.^(Zvector/10))./((distanceBS_To_DownlinkUsers_Zone1/RadiusOfNearestUser).^ploss);
    
    PL_downlink_Zone1 = GetPathloss(103.8, 20.9, distanceBS_To_DownlinkUsers_Zone1/1000);

        
    D_H1k = (diag(PL_downlink_Zone1)).^(0.5)
    
    
%% large-fading for L downlink users - outer zone

%     Zvector = StandardDeviation*randn(1,L);
    rvector = InnerZoneRadius*ones(1,L) + (RadiusOfCell-InnerZoneRadius)*rand(1,L);
    anglevector = 2*pi*rand(1,L);
    
    positionDownlinkUsers = [positionDownlinkUsers; (rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];

        
    distanceBS_To_DownlinkUsers_Zone2 = sqrt(sum(positionDownlinkUsers(K+1:K+L)'.^2));
        
%     betavector_downlink_Zone2 = (10.^(Zvector/10))./((distanceBS_To_DownlinkUsers_Zone2/RadiusOfNearestUser).^ploss);

    PL_downlink_Zone2 = GetPathloss(103.8, 20.9, distanceBS_To_DownlinkUsers_Zone2/1000);

        
    D_H2l = (diag(PL_downlink_Zone2)).^(0.5)
    
    
%% large-fading for P uplink users - outer zone

%     Zvector = StandardDeviation*randn(1,P);
    rvector = InnerZoneRadius*ones(1,P) + (RadiusOfCell-InnerZoneRadius)*rand(1,P);
    anglevector = 2*pi*rand(1,P);
    
    positionUplinkUsers = [(rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];

        
    distanceBS_To_UplinkUsers_Zone2 = sqrt(sum(positionUplinkUsers'.^2));
        
%     betavector_uplink_Zone2 = (10.^(Zvector/10))./((distanceBS_To_UplinkUsers_Zone2/RadiusOfNearestUser).^ploss);

    PL_uplink_Zone2 = GetPathloss(103.8, 20.9, distanceBS_To_UplinkUsers_Zone2/1000);

        
    D_G2p = (diag(PL_uplink_Zone2)).^(0.5)
    
    
%% large-fading for Q uplink users - inner zone

%     Zvector = StandardDeviation*randn(1,K);
    rvector = RadiusOfNearestUser*ones(1,Q) + (InnerZoneRadius-RadiusOfNearestUser)*rand(1,Q);
    anglevector = 2*pi*rand(1,Q);
    
    positionUplinkUsers = [positionUplinkUsers; (rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];
        
    distanceBS_To_UplinkUsers_Zone1 = sqrt(sum(positionUplinkUsers(P+1:P+Q)'.^2));
        
%     betavector_downlink_Zone1 = (10.^(Zvector/10))./((distanceBS_To_DownlinkUsers_Zone1/RadiusOfNearestUser).^ploss);
    
    PL_uplink_Zone1 = GetPathloss(103.8, 20.9, distanceBS_To_UplinkUsers_Zone1/1000);

        
    D_G1q = (diag(PL_uplink_Zone1)).^(0.5)
    
    
%% large-fading for CCI P ULUs to K DLUs 

    positionDownUsers_PQtimes = kron(positionDownlinkUsers(1:K+L,:),ones(P+Q,1));
    positionUpUsers_KLtimes = repmat(positionUplinkUsers(1:P+Q,:), K+L, 1);
    
    distance_PQup2KLdown = sqrt(sum((positionDownUsers_PQtimes-positionUpUsers_KLtimes).^2 ,2));
    
%     betavector_updownusers = (10.^(Zvector/10))./((distance_UpDownUsers/RadiusOfNearestUser).^ploss);

    PL_PQup2KLdown = GetPathloss(145.4, 37.5, distance_PQup2KLdown/1000);
    
    D_f2_pq_kl = (reshape(PL_PQup2KLdown, P+Q, K+L)).^(0.5)

%     positionDownUsers_Ptimes = kron(positionDownlinkUsers(1:K,:),ones(P,1));
%     positionUpUsers_Ktimes = repmat(positionUplinkUsers(1:P,:), K, 1);
%     
%     distance_Pup2Kdown = sqrt(sum((positionDownUsers_Ptimes-positionUpUsers_Ktimes).^2 ,2));
%     
% %     betavector_updownusers = (10.^(Zvector/10))./((distance_UpDownUsers/RadiusOfNearestUser).^ploss);
% 
%     PL_Pup2Kdown = GetPathloss(145.4, 37.5, distance_Pup2Kdown/1000);
%     
%     D_f2pk = (reshape(PL_Pup2Kdown, P, K)).^(0.5)
    
    
%% large-fading for CCI Q ULUs to L DLUs 

    D_f1ql = D_f2_pq_kl(P+1:P+Q,K+1:K+L)

%     positionDownUsers_Qtimes = kron(positionDownlinkUsers(K+1:K+L,:),ones(Q,1));
%     positionUpUsers_Ltimes = repmat(positionUplinkUsers(P+1:P+Q,:), L, 1);
%     
%     distance_Qup2Ldown = sqrt(sum((positionDownUsers_Qtimes-positionUpUsers_Ltimes).^2 ,2));
%     
% %     betavector_updownusers = (10.^(Zvector/10))./((distance_UpDownUsers/RadiusOfNearestUser).^ploss);
% 
%     PL_Qup2Ldown = GetPathloss(145.4, 37.5, distance_Qup2Ldown/1000);
%     
%     D_f1ql = (reshape(PL_Qup2Ldown, Q, L)).^(0.5)

      
    
%% large-fading for M Eavs - M1 Eavs in Zone 1

    rvector = RadiusOfNearestUser*ones(1,M1) + (InnerZoneRadius-RadiusOfNearestUser)*rand(1,M1);
    anglevector = 2*pi*rand(1,M1);
    
    positionEavs = [(rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];

    M2 = M - M1; % M2 Eavs in Zone 2
    
    rvector = InnerZoneRadius*ones(1,M2) + (RadiusOfCell-InnerZoneRadius)*rand(1,M2);
    anglevector = 2*pi*rand(1,M2);
    
    positionEavs = [positionEavs; (rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];
    
    % Pathloss BS to Eavs
    distanceBS_To_Eavs = sqrt(sum(positionEavs'.^2));
    
    D_Hm = (GetPathloss(103.8, 20.9, distanceBS_To_Eavs/1000)).^(0.5);
    
    
    % Pathlos ULUs to Eavs
    
    positionEavs_PplusQtimes = kron(positionEavs,ones(P+Q,1));
    positionUpUsers_Mtimes = repmat(positionUplinkUsers, M, 1);
    
    distance_UpUsers2Eavs = sqrt(sum((positionEavs_PplusQtimes-positionUpUsers_Mtimes).^2 ,2));
    
    PL_UpUsers2Eavs = GetPathloss(103.8, 20.9, distance_UpUsers2Eavs/1000);
    
    PL_UpUsers2Eavs_Mat = (reshape(PL_UpUsers2Eavs, P+Q, M)).^(0.5);
    
    D_Gm2p = PL_UpUsers2Eavs_Mat(1:P,:);
    
    D_Gm1q = PL_UpUsers2Eavs_Mat(P+1:P+Q,:);
    
    

end


function [Pathloss] = GetPathloss(firstPar, secondPar, distance)

% distance (km)

Pathloss_dB = firstPar + secondPar*log10(distance);

Pathloss = 10.^(-Pathloss_dB/10);

end

