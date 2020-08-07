close all
clear
clc

pathCell = regexp(path, pathsep, 'split');
if (~any(strcmpi('.\yalmip', pathCell)))
    disp('_________________ Adding path ... _________________')
    addpath(genpath('./yalmip'));
    addpath(genpath('./SeDuMi_1_3'));
    addpath(genpath('./SDPT3-4.0'));
end
clear pathCell

pathname = fileparts('./Figures/Rate_vs_EH/Layout_1Sum/');
addpath(genpath('./Figures/Rate_vs_EH/Layout_1Sum/'));

tic

%% Initialization

FileOrders = 31;

MethodUsed = 0; % 0 : Proposed method
                % 1 : Conven. FD method
                % 2 : Conven. FD with NOMA
                % 3 : HD
                % 4 : Proposed without Jamming

NumberOfRunning = 10;

ReusedChannel = 0;
Reusedfrom = 0;

layoutfile = ['layoutForAll_1.mat'];
loadingfile = [];%['[test] SRvsPbs_Conven.mat'];%['[test] SRvsPbs_HD.mat'];%['SuccessfullySolved1.mat'];

EavsOn = 1;

aux = [];
                
Methodname = {'Proposed'; 'ConvenWithoutNOMA'; 'ConvenWithNOMA'; 'HD'; 'ProposedNoJamming'};

if (MethodUsed==1 || MethodUsed==2)
    Conven = 1;
else
    Conven = 0;
end
Fixed_Timegroup = []; % an empty vector: for 2 time-groups optimization, 
                      % 1x2 vector     : fix time for 2 groups as
                      % corresponding to 2 values in the vector

if (Conven)
    Fixed_Timegroup = [1-10^(-2), 10^(-2)];
end

if (MethodUsed==3)
    Fixed_Timegroup = [1, 1];
end

sigma_SI_dB = -75;

MaxIteration = 50;

N_t = 5; % number of Tx antennas at BS
N_r = 5; % number of Rx antennas at BS

G = 2; % 2 regions - inner and outer zone

K = 2; % number of DLUs in inner zone
L = 2; % number of DLUs in outer zone

Q = 2; % number of ULUs in inner zone
P = 2; % number of ULUs in outer zone

M = 2; % number of Eavs

M1= 1; % number of Eavs in inner zone

N_Eavs = [2 2]; % number of antennas at Eavs -- length(N_Eavs) = M


if (length(N_Eavs) ~= M)
    disp('Error: length(N_Eavs) must equal M !!!');
    pause(10);
end

% Em_dBm = -18;

% Em = 10.^(Em_dBm/10);

% Em = [0:0.002:0.01 0.015 0.02:0.0025:0.03]; % EH per device
Em = [0:0.01:0.1]; % EH sum

% Em = [0:0.0002:0.001];

SumEH = 1;
delta = 0.5;

Pbs_dB = 46; % dowlink power [dBm]

Pbs = 10.^(Pbs_dB/10); % downlink power

P_2p_max = 23*ones(1,P); % uplink power zone 2 [dBm]
P_1q_max = 23*ones(1,Q); % uplink power zone 1 [dBm]

P_up_max_dB = [P_2p_max P_1q_max];
P_up_max = 10.^(P_up_max_dB/10);

if (MethodUsed==3)
    sigma_SI = zeros(size(sigma_SI_dB));
else
    sigma_SI = 10.^(sigma_SI_dB/10);
end


Noise = -174; % dBm/Hz
BW = 10; % MHz

sigma_noise = 10^(-5);%10^(Noise/10)*BW*10^6;%



% sigma_k = 0.01; %10^(Noise/10)*BW*10^6;
% 
% sigma_m = 0.01; %10^(Noise/10)*BW*10^6;

%% Set users' positions

if (isempty(layoutfile))

RadiusOfCell = 100;
InnerZoneRadius = 50;
RadiusOfNearestUser = 10;
StandardDeviation = 8;
ploss = 3;

Parameters = [RadiusOfCell InnerZoneRadius RadiusOfNearestUser StandardDeviation ploss];

[D_H1k, D_H2l, D_G2p, D_G1q, D_f2_pq_kl, D_F1ql, D_Hm, D_Gm2p, D_Gm1q, positionDownlinkUsers, positionUplinkUsers, positionEavs] = CreateLargeScaleFading( K, L, P, Q, M, M1, Parameters );

else
    
    load(layoutfile, 'D_H1k');
    load(layoutfile, 'D_H2l');
    load(layoutfile, 'D_G2p');
    load(layoutfile, 'D_G1q');
    load(layoutfile, 'D_f2_pq_kl');
    load(layoutfile, 'D_F1ql');
    load(layoutfile, 'D_Hm');
    load(layoutfile, 'D_Gm2p');
    load(layoutfile, 'D_Gm1q');
    load(layoutfile, 'positionDownlinkUsers');
    load(layoutfile, 'positionUplinkUsers');
    load(layoutfile, 'positionEavs');
    load(layoutfile, 'RadiusOfCell');
    load(layoutfile, 'InnerZoneRadius');

end

Layout = Plot_Layout( RadiusOfCell, InnerZoneRadius, positionDownlinkUsers, positionUplinkUsers, positionEavs );

%%
h = waitbar(0, 'Please wait ...','Name','Percent of completion');

for FileOrder = FileOrders
    

                
filename = ['[' Methodname{MethodUsed+1} aux '] Rate_vs_EH_' num2str(FileOrder) '.mat'];
savedfile = fullfile(pathname, filename);

if (ReusedChannel)
    loadingfile = ['[' Methodname{Reusedfrom+1} aux '] Rate_vs_EH_' num2str(FileOrder) '.mat'];
end



%% Create Channels

if (isempty(loadingfile))

    AllChannels = cell(NumberOfRunning, 9);
    D_Hm_mat = cell(NumberOfRunning, M);
    D_Gm2p_mat = cell(NumberOfRunning, M);
    D_Gm1q_mat = cell(NumberOfRunning, M);

    for i_NumOfSim=1:1:NumberOfRunning

        AllChannels{i_NumOfSim, 1} = CreateSmallScaleFading(1, 1, N_t+N_r, K)*D_H1k; % H1k
        AllChannels{i_NumOfSim, 2} = CreateSmallScaleFading(1, 1, N_t+N_r, L)*D_H2l; % H2l
        AllChannels{i_NumOfSim, 3} = CreateSmallScaleFading(1, 1, N_r+N_t, P)*D_G2p; % G2p
        AllChannels{i_NumOfSim, 4} = CreateSmallScaleFading(1, 1, N_r+N_t, Q)*D_G1q; % G1q
        AllChannels{i_NumOfSim, 5} = CreateSmallScaleFading(1, 1, P+Q, K+L).*D_f2_pq_kl; % D_f2_pq_kl F2pk
        AllChannels{i_NumOfSim, 6} = AllChannels{i_NumOfSim, 5}(P+1:P+Q,K+1:K+L); % F1ql
%         AllChannels{i_NumOfSim, 5} = CreateSmallScaleFading(1, 1, P, K).*D_F2pk; % D_f2_pq_kl F2pk
%         AllChannels{i_NumOfSim, 6} = CreateSmallScaleFading(1, 1, Q, L).*D_F1ql; % F1ql

        for iM = 1:1:M
            D_Hm_mat{i_NumOfSim, iM} = repmat(D_Hm(iM), N_t+N_r, N_Eavs(iM));
            AllChannels{i_NumOfSim, 7}{iM} = CreateSmallScaleFading(1, 1, N_t+N_r, N_Eavs(iM)).*D_Hm_mat{i_NumOfSim, iM}; % Hm

            D_Gm2p_mat{i_NumOfSim, iM} = diag(D_Gm2p(:,iM));
            AllChannels{i_NumOfSim, 8}{iM} = D_Gm2p_mat{i_NumOfSim, iM}*CreateSmallScaleFading(1, 1, P, N_Eavs(iM)); % Gm2p

            D_Gm1q_mat{i_NumOfSim, iM} = diag(D_Gm1q(:,iM));
            AllChannels{i_NumOfSim, 9}{iM} = D_Gm1q_mat{i_NumOfSim, iM}*CreateSmallScaleFading(1, 1, Q, N_Eavs(iM)); % Gm1q
        end

        AllChannels{i_NumOfSim, 10} = CreateSmallScaleFading(1, 1, N_t, N_r); % G_SI

    end

else
    
    load(loadingfile,'AllChannels');
    load(loadingfile,'RadiusOfCell');
    load(loadingfile,'InnerZoneRadius');
    load(loadingfile,'positionDownlinkUsers');
    load(loadingfile,'positionUplinkUsers');
    load(loadingfile,'positionEavs');

end


XDir = Em;

Alg1_Rate_all = zeros(NumberOfRunning, length(XDir));

DownlinkRate_PerGroupPerUser = zeros(K, G, NumberOfRunning);


% Clus = parcluster('local');
% Clus.NumWorkers = 4;
% 
% poolobj = parpool(Clus, Clus.NumWorkers);

% i_NumOfSim = 0;

for i_NumOfSim=1:1:NumberOfRunning
    
    % assign channels
    
    if (MethodUsed~=3)
        H1k = AllChannels{i_NumOfSim, 1}(1:N_t,:);
        H2l = AllChannels{i_NumOfSim, 2}(1:N_t,:);
        G2p = AllChannels{i_NumOfSim, 3}(1:N_r,:);
        G1q = AllChannels{i_NumOfSim, 4}(1:N_r,:);
    else
        H1k = AllChannels{i_NumOfSim, 1};
        H2l = AllChannels{i_NumOfSim, 2};
        G2p = AllChannels{i_NumOfSim, 3};
        G1q = AllChannels{i_NumOfSim, 4};
    end
    if (Conven)
        F2pk = AllChannels{i_NumOfSim, 5};
        F1ql = [];
    else
        if (MethodUsed==3)
            F2pk = zeros(K+L,P+Q);
            F1ql = zeros(K+L,P+Q);
        else
            F2pk = AllChannels{i_NumOfSim, 5}(1:P,1:K);
            F1ql = AllChannels{i_NumOfSim, 6};
        end
    end
    if (EavsOn)
        Hm = cell(1,M);
        Gm2p = cell(1,M);
        Gm1q = cell(1,M);
        for iM = 1:1:M
            if (MethodUsed~=3)
                Hm{iM} = AllChannels{i_NumOfSim, 7}{iM}(1:N_t,:);
            else
                Hm{iM} = AllChannels{i_NumOfSim, 7}{iM};
            end
            
            if (Conven)
                Gm2pq{iM} = [AllChannels{i_NumOfSim, 8}{iM}; AllChannels{i_NumOfSim, 9}{iM}];
                Gm1q{iM} = [];
            else
                if (MethodUsed==3)
                    Gm2p{iM} = [];
                    Gm1q{iM} = [AllChannels{i_NumOfSim, 8}{iM}; AllChannels{i_NumOfSim, 9}{iM}];
                else
                    Gm2p{iM} = AllChannels{i_NumOfSim, 8}{iM};
                    Gm1q{iM} = AllChannels{i_NumOfSim, 9}{iM};
                end
            end

        end
    else
        Hm = cell(1,0);
        Gm2p = cell(1,0);
        Gm1q = cell(1,0);
    end
    
    if (MethodUsed~=3)
        G_SI = AllChannels{i_NumOfSim, 10}; 
    else
        G_SI = zeros(N_t+N_r,N_t+N_r);
    end
    
    if (Conven)
        
        AllChannel = {[H1k, H2l], [], [G2p, G1q], [], F2pk, F1ql, Hm, Gm2pq, Gm1q, G_SI};
        
    else
        if (MethodUsed==3)
            AllChannel = {[H1k, H2l], [], [], [G2p, G1q], F2pk, F1ql, Hm, Gm2p, Gm1q, G_SI};
        else
            AllChannel = {H1k, H2l, G2p, G1q, F2pk, F1ql, Hm, Gm2p, Gm1q, G_SI};
        end
    
    end
    % notify completed percent
    
    percent = ((FileOrder-FileOrders(1))*NumberOfRunning + (i_NumOfSim-1)) / (length(FileOrders)*NumberOfRunning);
    waitbar(percent, h,[ num2str(floor(percent*100)) ' % Completed'])

Clus = parcluster('local');
Clus.NumWorkers = 6;

poolobj = parpool(Clus, Clus.NumWorkers);
    
Alg1_OptValue = zeros(1, length(XDir));
EH = zeros(length(XDir), M+1);

% for iXDir = 1:1:length(XDir)
parfor iXDir = 1:1:length(XDir)


    
    %% Algorithm 1
    if (SumEH)
        EH(iXDir,:) = [SumEH delta XDir(iXDir)];
    else
        EH(iXDir,:) = [SumEH delta XDir(iXDir)*ones(1,M)];
    end
    [isBreak1, Alg1_OptValue(iXDir), Alg1_OptValueChain, OptSolution ] = ProposedAlg1( Pbs, P_up_max, AllChannel, [sigma_SI sigma_noise], MaxIteration, Fixed_Timegroup, MethodUsed, 0, [], EH(iXDir,:)); % P_bs = XDir(iXDir)
%     Alg1_OptValue = 0;
%     
%     if (isBreak1)
%         NumOfSamples1 = NumOfSamples1 - 1;
%         i_NumOfSim = i_NumOfSim - 1;
%         break;
%     end
    
    
    Alg1_Rate_all(i_NumOfSim, iXDir) = Alg1_OptValue(iXDir);

    
    
    
end

%     waitbar(i_NumOfSim / NumberOfRunning, h,[ str_Group ':  ' num2str(floor(i_NumOfSim*100 / NumberOfRunning)) ' % Completed'])
    
delete(poolobj);

clear Clus


end

% delete(poolobj);
% 
% clear Clus



%%

if (NumberOfRunning == 1)
    Alg1_Rate = Alg1_Rate_all
%     Alg2_Rate = Alg2_Rate_all
%     Alg3_Rate = Alg3_Rate_all
else
    Alg1_Rate = mean(Alg1_Rate_all)
%     Alg1_Rate = mean(Alg1_Rate_all(1:5,:))
%     Alg2_Rate = mean(Alg2_Rate_all)
%     Alg3_Rate = mean(Alg3_Rate_all)
end

% HD_Rate = repmat(mean(OptValue_HalfDuplex_Stat), 1, length(rho_dB));

save(savedfile,'-regexp','^(?!h$).');


end

close(h);

%% Plot

figure;
% openfig('[test] SRvsPbs_Conven.fig');

% hold on
% 
plot(Em, Alg1_Rate, 'k-', 'linewidth', 2, 'markersize',9);


% if (XDir==1)
%     plot([1:1:length(Alg1_OptValueChain)], Alg1_OptValueChain, 'b-', 'linewidth', 2, 'markersize',9);
% else
%     plot(XDir, Alg1_Rate, 'bs--', 'linewidth', 2, 'markersize',9);
% end
% plot([1:1:length(Alg2_OptValueChain_sum)], Alg2_OptValueChain_sum, 'r-', 'linewidth', 2, 'markersize',9);
% plot([1:1:length(Alg3_OptValueChain)], Alg3_OptValueChain, 'g-', 'linewidth', 2, 'markersize',9);
% plot([1:1:length(Alg3_OptValueChain)], Alg3_OptValueChain*eta, 'g-', 'linewidth', 2, 'markersize',9);


% plot(full_range_rho, Alg1_Rate, 'bs--', 'linewidth', 2, 'markersize',9);
% plot(full_range_rho, Alg2_Rate, 'rs-', 'linewidth', 2, 'markersize',9);
% plot(rho_dB, HD_Rate, 'k^-', 'linewidth', 2, 'markersize',9);

% legend('FD - fixed time group', 'FD - optimal time group');

time = toc
hour = 0; min=0;
if (time>3600)
    hour = floor(time/3600);
    time = mod(time,3600);
end
if (time>60)
    min = floor(time/60);
    time = mod(time,60);
end
disp(['Running Time = ' num2str(hour) ' h ' num2str(min) ' m ' num2str(time) ' s.' ]);