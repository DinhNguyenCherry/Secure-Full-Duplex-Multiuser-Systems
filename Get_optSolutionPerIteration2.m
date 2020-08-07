function [ OptimalValue, opt_X, opt_alpha, opt_beta, opt_mu, alpha_p_next, beta_p_next, beta_SCSI_next, Status ] = Get_optSolutionPerIteration2(Pbs, P_up_max, AllChannels, sigma, Fixed_Timegroup, X_current, alpha_current, beta_current, mu_current, alpha_p_current, beta_p_current, preStatus, Method, Ru, beta_SCSI_current, eps_SCSI, EH)
%GET_OPTSOLUTIONPERITERATION Summary of this function goes here
%   Detailed explanation goes here

%% capture parameters

if (nargin < nargin('Get_optSolutionPerIteration2'))
    beta_SCSI_current =  {[], [], [], []};
    eps_SCSI = [];
end



K = size(AllChannels{1, 1}, 2);
L = size(AllChannels{1, 2}, 2);
P = size(AllChannels{1, 3}, 2);
Q = size(AllChannels{1, 4}, 2);
M = size(AllChannels{1, 7}, 2);

N_t = size(AllChannels{1, 1}, 1);
N_r = size(AllChannels{1, 3}, 1);

if (~isempty(EH))
    SumEH = EH(1);
    delta = EH(2);
    if (SumEH)
        Em = EH(3);
    else
        Em = EH(3:end);
    end
else
    SumEH = 0;
end

if (length(Method)>1)
    IsDownlink = Method(2);
    Method = 3;
    if (IsDownlink)
        L = 0;
        P = 0;
        Q = 0;
    else
        K = 0;
        L = 0;
        P = 0;
    end
end

if (length(eps_SCSI)>1)
    WCS = 0;
    SCSI = 1;
    
    beta_1k_DL_current = beta_SCSI_current{1,1};
    beta_2l_DL_current = beta_SCSI_current{1,2};
    beta_1q_UL_current = beta_SCSI_current{1,3};
    beta_2p_UL_current = beta_SCSI_current{1,4};
    
    eps_m1k = eps_SCSI(1);
    eps_m2l = eps_SCSI(2);
    eps_m1q = eps_SCSI(3);
    eps_m2p = eps_SCSI(4);
else
    if (length(eps_SCSI)==1)
        WCS = 1;
        SCSI = 0;
    else
        WCS = 0;
        SCSI = 0;
    end
    
end

%% assign channels, SI and AWGN
    
H1k = AllChannels{1, 1};
H2l = AllChannels{1, 2};
G2p = AllChannels{1, 3};
G1q = AllChannels{1, 4};
F2pk = AllChannels{1, 5};
F1ql = AllChannels{1, 6};
Hm = cell(1,M);
Gm2p = cell(1,M);
Gm1q = cell(1,M);
for iM = 1:1:M

    Hm{iM} = AllChannels{1, 7}{iM};

    Gm2p{iM} = AllChannels{1, 8}{iM};

    Gm1q{iM} = AllChannels{1, 9}{iM};

end

G_SI = AllChannels{1, 10};


sigma_SI = sigma(1);
sigma_noise = sigma(2); %0.01*ones(1, K);


%% extract current solution

X1 = X_current{1,1};
X2 = X_current{1,2};

W1_current = X1{1,1};
V1_current = X1{1,2};
rho2p_current = X1{1,3};

W2_current = X2{1,1};
V2_current = X2{1,2};
rho1q_current = X2{1,3};

mu_m1k_current = mu_current{1,1};
mu_m2l_current = mu_current{1,2};
mu_tilde_m1q_current = mu_current{1,3};
mu_tilde_m2p_current = mu_current{1,4};


%% initialize variables

W1_next = sdpvar(N_t, K, 'full','complex');
if (Method==4)
    V1_next = zeros(N_t,N_t);
else
if (Method~=3)
    V1_next = sdpvar(N_t, N_t, 'full','complex');%zeros(N_t,N_t);;%sdpvar(N_t, N_t, 'full','complex');%
else
    if (IsDownlink)
        V1_next = sdpvar(N_t, N_t, 'full','complex');%
    else
        V1_next = zeros(N_t,N_t);%
    end
end
end
rho2p_next = sdpvar(1, P, 'full','real');

W2_next = sdpvar(N_t, L, 'full','complex');
% V2_next = sdpvar(N_t, N_t, 'full','complex');
if (Method) % Method==3||Method==1||Method==2
    V2_next = zeros(N_t,N_t);
else
    V2_next = sdpvar(N_t, N_t, 'full','complex');%zeros(N_t,N_t);%
end
rho1q_next = sdpvar(1, Q, 'full','real');

if (Method==4)
    Method=0;
end

eta = sdpvar(1, 1, 'full','real');

GammaDL_1k = sdpvar(1, K, 'full','real');
GammaDL_2l = sdpvar(1, L, 'full','real');
GammaUL_1q = sdpvar(1, Q, 'full','real');
GammaUL_2p = sdpvar(1, P, 'full','real');

% GammaDL_1k = zeros(1,K);
% GammaDL_2l = zeros(1,L);
% GammaUL_1q = zeros(1,Q);
% GammaUL_2p = zeros(1,P);

if (~isempty(Fixed_Timegroup))
    alpha_next = 1/Fixed_Timegroup(1);
    beta_next = 1/Fixed_Timegroup(2);
else
    alpha_next = sdpvar(1, 1, 'full','real');
    beta_next = sdpvar(1, 1, 'full','real');
end

mu_m1k_next = sdpvar(M, K, 'full','real');
mu_m2l_next = sdpvar(M, L, 'full','real');
mu_tilde_m1q_next = sdpvar(M, Q, 'full','real');
mu_tilde_m2p_next = sdpvar(M, P, 'full','real');

% mu_m1k_next = zeros(M, K);
% mu_m2l_next = zeros(M, L);
% mu_tilde_m1q_next = zeros(M, Q);
% mu_tilde_m2p_next = zeros(M, P);

%% add objective function and constraints

obj = eta;
cons = [];

% cons = [cons, eta>= 10^(-5)];

cons = [cons, GammaDL_1k>=10^(-30)];
cons = [cons, GammaUL_1q>=10^(-30)];
if (Method~=3)
cons = [cons, GammaDL_2l>=10^(-30)];

cons = [cons, GammaUL_2p>=10^(-30)];
end


%% Constraints (13e)

for iQ = 1:1:Q
    cons = [cons, rho1q_next(iQ) >= 10^(-30)];
end

for iP = 1:1:P
    cons = [cons, rho2p_next(iP) >= 10^(-30)];
end


%% Constraints (17)


if (isempty(Fixed_Timegroup))

    cons = [cons, alpha_next>=1+10^(-4)];
    cons = [cons, beta_next>=1+10^(-4)];
    cons = [cons, alpha_next<=100];
    cons = [cons, beta_next<=100];
    
    alpha_p = sdpvar(1, 1, 'full','real');
    beta_p = sdpvar(1, 1, 'full','real');

    cons = [cons, cone([1, 0.5*(alpha_next - alpha_p)], 0.5*(alpha_next + alpha_p))];
    cons = [cons, cone([1, 0.5*(beta_next - beta_p)], 0.5*(beta_next + beta_p))];
    cons = [cons, alpha_p>=10^(-2)];
    cons = [cons, beta_p>=10^(-2)];

    % tau = sdpvar(1,1);
    % cons = [cons, tau>=10^(-4)];
    % cons = [cons, tau<=1];
    % cons = [cons, cone([1 0.5*(tau-alpha_next)], 0.5*(tau+alpha_next))];
    % cons = [cons, cone([1 0.5*(1-tau-beta_next)], 0.5*(1-tau+alpha_next))];

    cons = [cons, (alpha_p + beta_p) <= 1 ];
    % cons = [cons, (alpha_p + beta_p) >= 1 ];
else
    
    alpha_p = Fixed_Timegroup(1);
    beta_p = Fixed_Timegroup(2);

end





%% Constraints (20d)
if (Method==3)
    
    for iQ = 1:1:Q    
        cons = [cons, cone( rho1q_next(iQ), sqrt(P_up_max(P+iQ)) ) ];   
    end
else
    for iQ = 1:1:Q    
        cons = [cons, cone([rho1q_next(iQ), 0.5*(beta_next-P_up_max(P+iQ))], 0.5*(beta_next+P_up_max(P+iQ))) ];   
    end
end

%% Constraints (25)

for iK = 1:1:K
    cons = [cons, real(H1k(:,iK)'*W1_current(:,iK))*(2*real(H1k(:,iK)'*W1_next(:,iK)) - real(H1k(:,iK)'*W1_current(:,iK))) >= 0 ];
end


%% Constraints (27)

A1k_mathcal = zeros(1,K);
B1k_mathcal = zeros(1,K);
C1k_mathcal = zeros(1,K);

varphi1k_p = sdpvar(1, K, 'full','real');
varphi1k_bar = sdpvar(1, K, 'full','real');
C1k_DL = sdpvar(1, K, 'full','real');
if (Method==2)
    K2 = floor(K/2);
    A1k_mathcal_p = zeros(1,K2);
    B1k_mathcal_p = zeros(1,K2);
    C1k_mathcal_p = zeros(1,K2);
    
    varphi1k_pp = sdpvar(1, K2, 'full','real');
    varphi1k_bar_p = sdpvar(1, K2, 'full','real');
    C1k_DL_p = sdpvar(1, K, 'full','real');
end



for iK = 1:1:K
    if (Method==2 && iK>K2)
%         if (abs(H1k(:,iK)'*W1_current(:,iK))^2<=abs(H1k(:,iK-K2)'*W1_current(:,iK))^2)
%             h_channel_min = H1k(:,iK);
%         else
%             h_channel_min = H1k(:,iK-K2);
%         end
        [A1k_mathcal(iK), B1k_mathcal(iK), C1k_mathcal(iK)] = Get_ABCmathcal( X1, alpha_current, iK, H1k(:,iK), F2pk(:,iK), sigma_noise, Method );
        [A1k_mathcal_p(iK-K2), B1k_mathcal_p(iK-K2), C1k_mathcal_p(iK-K2)] = Get_ABCmathcal( X1, alpha_current, iK, H1k(:,iK-K2), F2pk(:,iK), sigma_noise, Method );
    else
        [A1k_mathcal(iK), B1k_mathcal(iK), C1k_mathcal(iK)] = Get_ABCmathcal( X1, alpha_current, iK, H1k(:,iK), F2pk(:,iK), sigma_noise, Method );
    end
    
    varphi1k_p(iK) = real(H1k(:,iK)'*W1_current(:,iK))*(2*real(H1k(:,iK)'*W1_next(:,iK)) - real(H1k(:,iK)'*W1_current(:,iK)));
    
    [~, varphi_vec] = Get_varphi( {W1_next, V1_next, rho2p_next}, iK, H1k(:,iK), F2pk(:,iK), sigma_noise, Method );
    
    cons = [cons, cone([varphi_vec 0.5*(varphi1k_p(iK) - varphi1k_bar(iK))], 0.5*(varphi1k_p(iK) + varphi1k_bar(iK))) ];
    
    if (Method==2 && iK>K2)
        varphi1k_pp(iK-K2) = 2*real(W1_current(:,iK)'*H1k(:,iK-K2)*H1k(:,iK-K2)'*W1_next(:,iK)) - abs(H1k(:,iK-K2)'*W1_current(:,iK))^2;
        [~, varphi_vec_p] = Get_varphi( {W1_next, V1_next, rho2p_next}, iK, H1k(:,iK-K2), F2pk(:,iK), sigma_noise, Method );
        cons = [cons, cone([varphi_vec_p 0.5*(varphi1k_pp(iK-K2) - varphi1k_bar_p(iK-K2))], 0.5*(varphi1k_pp(iK-K2) + varphi1k_bar_p(iK-K2))) ];
    end
    
    C1k_DL(iK) = real(A1k_mathcal(iK) - B1k_mathcal(iK)*varphi1k_bar(iK) - C1k_mathcal(iK)*alpha_next);
    
    cons = [cons, C1k_DL(iK)>=10^(-20)];
    
    cons = [cons, C1k_DL(iK) >= real(eta + GammaDL_1k(iK)) ];
    
    if (Method==2 && iK>K2)
        C1k_DL_p(iK-K2) = real(A1k_mathcal_p(iK-K2) - B1k_mathcal_p(iK-K2)*varphi1k_bar_p(iK-K2) - C1k_mathcal_p(iK-K2)*alpha_next);
    
        cons = [cons, C1k_DL_p(iK-K2)>=10^(-20)];

        cons = [cons, C1k_DL_p(iK-K2) >= real(eta + GammaDL_1k(iK)) ];
    end
    
%     cons = [cons, varphi1k_bar(iK)>=10^(-5)];
    
%     if (Ru)
%         if (M>0 && isempty(findstr(preStatus{3},'NoEav')))
%             cons = [cons, real(A1k_mathcal(iK) - B1k_mathcal(iK)*varphi1k_bar(iK) - C1k_mathcal(iK)*alpha_next) >= real(eta + GammaDL_1k(iK)) ];
%         else
%             cons = [cons, real(A1k_mathcal(iK) - B1k_mathcal(iK)*varphi1k_bar(iK) - C1k_mathcal(iK)*alpha_next) >= 0.1 ];
%         end
%     else
        
%     end
    
end


%% Constraints (28)

A2l_mathcal = zeros(1,L);
B2l_mathcal = zeros(1,L);
C2l_mathcal = zeros(1,L);
varphi2l_p = sdpvar(1, L, 'full','real');
varphi2l_bar = sdpvar(1, L, 'full','real');

C2l_DL = sdpvar(1, L, 'full','real');

for iL = 1:1:L
    
    [A2l_mathcal(iL), B2l_mathcal(iL), C2l_mathcal(iL)] = Get_ABCmathcal( X2, beta_current, iL, H2l(:,iL), F1ql(:,iL), sigma_noise, Method );
    
    varphi2l_p(iL) = real(H2l(:,iL)'*W2_current(:,iL))*(2*real(H2l(:,iL)'*W2_next(:,iL)) - real(H2l(:,iL)'*W2_current(:,iL)));
    
    [~, varphi_vec] = Get_varphi( {W2_next, V2_next, rho1q_next}, iL, H2l(:,iL), F1ql(:,iL), sigma_noise, Method );
    
    cons = [cons, cone([varphi_vec 0.5*(varphi2l_p(iL) - varphi2l_bar(iL))], 0.5*(varphi2l_p(iL) + varphi2l_bar(iL))) ];
    
    C2l_DL(iL) = A2l_mathcal(iL) - B2l_mathcal(iL)*varphi2l_bar(iL) - C2l_mathcal(iL)*beta_next;
    
    cons = [cons, C2l_DL(iL)>=10^(-20)];
    
%     if (Ru)
%         if (M>0 && isempty(findstr(preStatus{3},'NoEav')))
%             cons = [cons, real((A2l_mathcal(iL) - B2l_mathcal(iL)*varphi2l_bar(iL) - C2l_mathcal(iL)*beta_next)) >= real(eta + GammaDL_2l(iL)) ];
%         else
%             cons = [cons, real((A2l_mathcal(iL) - B2l_mathcal(iL)*varphi2l_bar(iL) - C2l_mathcal(iL)*beta_next)) >= 0.1 ];
%         end
%     else
        cons = [cons, real((A2l_mathcal(iL) - B2l_mathcal(iL)*varphi2l_bar(iL) - C2l_mathcal(iL)*beta_next)) >= real(eta + GammaDL_2l(iL)) ];
%     end
    
end


%% Constraints (30)
if (Method~=3)
    for iL = 1:1:L
        cons = [cons, real(H2l(:,iL)'*W2_current(:,iL))*(2*real(H2l(:,iL)'*W2_next(:,iL)) - real(H2l(:,iL)'*W2_current(:,iL))) >= 0 ];
    end
end


%% Constraints (33)

A1q_tilde_mathcal = zeros(1,Q);
B1q_tilde_mathcal = zeros(1,Q);
C1q_tilde_mathcal = zeros(1,Q);
phi1q_bar = sdpvar(1, Q, 'full','real');
phi1q_vec = cell(1,Q);

C1q_UL = sdpvar(1, Q, 'full','real');

for iQ = 1:1:Q
    [A1q_tilde_mathcal(iQ), B1q_tilde_mathcal(iQ), C1q_tilde_mathcal(iQ), Omega1q] = Get_ABCtildemathcal( X2, beta_current, iQ, G1q(:,iQ:end), G_SI, sigma_SI, sigma_noise );
    
    [phi1q_vec{iQ}] = Get_phivec( {W2_next, V2_next, rho1q_next}, iQ, G1q(:,iQ:end), G_SI, sigma_SI, sigma_noise, Omega1q );
    
    cons = [cons, cone([phi1q_vec{iQ} 0.5*(phi1q_bar(iQ) - beta_current)], 0.5*(phi1q_bar(iQ) + beta_current) ) ];
    
    C1q_UL(iQ) = A1q_tilde_mathcal(iQ) + B1q_tilde_mathcal(iQ)*rho1q_next(iQ) - phi1q_bar(iQ) - C1q_tilde_mathcal(iQ)*beta_next;
    
    
    
    if (Ru)

         cons = [cons, real(A1q_tilde_mathcal(iQ) + B1q_tilde_mathcal(iQ)*rho1q_next(iQ) - phi1q_bar(iQ) - C1q_tilde_mathcal(iQ)*beta_next)>= Ru + GammaUL_1q(iQ)];
    
    else
        
        cons = [cons, C1q_UL(iQ)>=10^(-20)];
        cons = [cons, real(A1q_tilde_mathcal(iQ) + B1q_tilde_mathcal(iQ)*rho1q_next(iQ) - phi1q_bar(iQ) - C1q_tilde_mathcal(iQ)*beta_next) >= real(eta + GammaUL_1q(iQ)) ];
        
    end
    
end


%% Constraints (34)

A2p_tilde_mathcal = zeros(1,P);
B2p_tilde_mathcal = zeros(1,P);
C2p_tilde_mathcal = zeros(1,P);
phi2p_bar = sdpvar(1, P, 'full','real');
phi2p_vec = cell(1,P);

C2p_UL = sdpvar(1, P, 'full','real');

for iP = 1:1:P
    [A2p_tilde_mathcal(iP), B2p_tilde_mathcal(iP), C2p_tilde_mathcal(iP), Omega2p] = Get_ABCtildemathcal( X1, alpha_current, iP, G2p(:,iP:end), G_SI, sigma_SI, sigma_noise );
    
    [phi2p_vec{iP}] = Get_phivec( {W1_next, V1_next, rho2p_next}, iP, G2p(:,iP:end), G_SI, sigma_SI, sigma_noise, Omega2p );
    
    cons = [cons, cone([phi2p_vec{iP} 0.5*(phi2p_bar(iP) - alpha_current)], 0.5*(phi2p_bar(iP) + alpha_current) ) ];
    
    C2p_UL(iP) = A2p_tilde_mathcal(iP) + B2p_tilde_mathcal(iP)*rho2p_next(iP) - phi2p_bar(iP) - C2p_tilde_mathcal(iP)*alpha_next;
    
    
    
    if (Ru)
        
        cons = [cons, real(A2p_tilde_mathcal(iP) + B2p_tilde_mathcal(iP)*rho2p_next(iP) - phi2p_bar(iP) - C2p_tilde_mathcal(iP)*alpha_next)>= Ru + GammaUL_2p(iP)];

    else
        
        cons = [cons, C2p_UL(iP)>=10^(-20)];
        cons = [cons, real(A2p_tilde_mathcal(iP) + B2p_tilde_mathcal(iP)*rho2p_next(iP) - phi2p_bar(iP) - C2p_tilde_mathcal(iP)*alpha_next) >= real(eta + GammaUL_2p(iP)) ];
    
    end
    
end



%% Constraints (37a) (39) (41) (42) for Eavs

% beta_dot = sdpvar(1,1,'full','real');
% cons = [cons, cone([1, 0.5*(beta_next-1 - beta_dot)], 0.5*(beta_next-1 + beta_dot)) ];

% if (~isempty(EH))
%     etaEH = sdpvar(1, 1, 'full','real'); 
%     amEH = sdpvar(1, M, 'full','real');
%     bmEH = sdpvar(1, M, 'full','real');
%     cmEH = zeros(1, M);
%     dmEH = zeros(1, M);
%     Lambda_approxEH = sdpvar(1,M, 'full', 'real');
%     EHm = sdpvar(1, M, 'full','real');
%     for iM = 1:1:M
%         [amEH(iM), bmEH(iM), cmEH(iM), dmEH(iM)] = Get_abcd_m({W1_next, V1_next, rho2p_next}, {W2_next, V2_next, rho1q_next}, X1, X2, Hm{iM}, Gm2p{iM}, Gm1q{iM}, SCSI); % for f_mathcal and g_mathcal
%     
%         Lambda_approxEH(iM) = (2/alpha_current)*amEH(iM) + (2/beta_current)*bmEH(iM) - cmEH(iM)*alpha_next/(alpha_current^2)-dmEH(iM)*beta_next/(beta_current^2);
% 
%  
%         disp('Add EH constraint hereeeee');
%         
%         EHm(iM) = delta*Lambda_approxEH(iM);
% %         cons = [cons, EHm(iM)>=Em(iM)];
% %         cons = [cons, EHm(iM)-Em(iM)>=etaEH];
%         cons = [cons, EHm(iM)>=etaEH];
% 
%     end
% end

if (M>0 && isempty(findstr(preStatus{3},'NoEav')))
    
    disp('Add the constraints associated to EAVS hereeeeeeeee');
    
am = sdpvar(1, M, 'full','real');
bm = sdpvar(1, M, 'full','real');
cm = zeros(1, M);
dm = zeros(1, M);

mu_p_m1k = sdpvar(M, K, 'full','real');
mu_p_m2l = sdpvar(M, L, 'full','real');
mu_tilde_p_m1q = sdpvar(M, Q, 'full','real');
mu_tilde_p_m2p = sdpvar(M, P, 'full','real');

mu_hat_m1k = sdpvar(M, K, 'full','real');
mu_hat_m2l = sdpvar(M, L, 'full','real');
mu_hat_m1q = sdpvar(M, Q, 'full','real');
mu_hat_m2p = sdpvar(M, P, 'full','real');

Lambda1 = zeros(1, M);
Lambda2 = zeros(1, M);

if (WCS)
    gammahat_ED_1k = zeros(M, K);
    gammahat_ED_2l = zeros(M, L);
    gammahat_EU_1q = zeros(M, Q);
    gammahat_EU_2p = zeros(M, P);
    
%     Theta_ED_1k = cell(M, K);
%     Theta_ED_2l = cell(M, L);
%     Theta_tilde_EU_1q = cell(M, Q);
%     Theta_tilde_EU_2p = cell(M, P);
    
    Theta_ED_1k = sdpvar(M, K, 'full','real');
    Theta_ED_2l = sdpvar(M, L, 'full','real');
    Theta_tilde_EU_1q = sdpvar(M, Q, 'full','real');
    Theta_tilde_EU_2p = sdpvar(M, P, 'full','real');
    
    C_hat_EDL_1k = sdpvar(M, K, 'full','real');
    C_hat_EDL_2l = sdpvar(M, L, 'full','real');
    C_hat_EUL_1q = sdpvar(M, Q, 'full','real');
    C_hat_EUL_2p = sdpvar(M, P, 'full','real');
    
else
    gamma1k_EDL = zeros(M, K);
    gamma2l_EDL = zeros(M, L);
    gamma1q_EUL = zeros(M, Q);
    gamma2p_EUL = zeros(M, P);

    C_tilde_EDL_1k = sdpvar(M, K, 'full','real');
    C_tilde_EDL_2l = sdpvar(M, L, 'full','real');
    C_tilde_EUL_1q = sdpvar(M, Q, 'full','real');
    C_tilde_EUL_2p = sdpvar(M, P, 'full','real');
    
    if (~isempty(EH))
        EHm = sdpvar(1, M, 'full','real');
    end
end

F_mathcal_1k = sdpvar(M, K, 'full', 'real');
G_mathcal_1q = sdpvar(M, Q, 'full', 'real');
G_mathcal_2l = sdpvar(M, L, 'full', 'real');
F_mathcal_2p = sdpvar(M, P, 'full', 'real');

N_em = zeros(1,M);
Lambda_approx = sdpvar(1,M, 'full', 'real');

if (WCS)
    XI = cell(2,M);
    if ((Method~=3) || (Method==3 && IsDownlink))
    SumJamming = cell(1,M);
    SumJamming_approx = sdpvar(N_t,N_t, 'hermitian', 'complex');
    SumJamming_approx_alpha = sdpvar(N_t,N_t, 'hermitian', 'complex');
    SumJamming_approx_beta = sdpvar(N_t,N_t, 'hermitian', 'complex');

    SumJamming_approx_alpha = (1/alpha_current)*(V1_next*V1_current' + V1_current*V1_next' - (1/alpha_current)*alpha_next*(V1_current*V1_current'));
    SumJamming_approx_beta = (1/beta_current)*(V2_next*V2_current' + V2_current*V2_next' - (1/beta_current)*beta_next*(V2_current*V2_current'));
    SumJamming_approx = (1/alpha_current)*(V1_next*V1_current' + V1_current*V1_next' - (1/alpha_current)*alpha_next*(V1_current*V1_current')) ...
                            + (1/beta_current)*(V2_next*V2_current' + V2_current*V2_next' - (1/beta_current)*beta_next*(V2_current*V2_current'));
    else
        SumJamming_approx_alpha = 0;
        SumJamming_approx_beta = 0;
        SumJamming_approx = 0;
    end

else

if (SCSI)
    beta_1k_DL = sdpvar(1, K, 'full', 'real');
    beta_2l_DL = sdpvar(1, L, 'full', 'real');
    beta_1q_UL = sdpvar(1, Q, 'full', 'real');
    beta_2p_UL = sdpvar(1, P, 'full', 'real');
    
    beta_bar_1k_DL = sdpvar(M, K, 'full', 'real');
    beta_bar_2l_DL = sdpvar(M, L, 'full', 'real');
    beta_bar_1q_UL = sdpvar(M, Q, 'full', 'real');
    beta_bar_2p_UL = sdpvar(M, P, 'full', 'real');
    
    xi_beta_1k_DL  = sdpvar(1, K, 'full', 'real');
    xi_beta_2l_DL  = sdpvar(1, L, 'full', 'real');
    xi_beta_1q_UL  = sdpvar(1, Q, 'full', 'real');
    xi_beta_2p_UL  = sdpvar(1, P, 'full', 'real');
end
end

for iM = 1:1:M
    
    N_em(iM) = size(Hm{iM},2);
    
    if (WCS)
        if ((Method~=3) || (Method==3 && IsDownlink))
            SumJamming{iM} = Get_SumJamming(V1_current, V2_current, alpha_current, beta_current, Hm{iM});
            XI{1,iM} = alpha_current*(SumJamming{iM} + sigma_noise^2*eye(N_em(iM)));
            XI{2,iM} = beta_current*(SumJamming{iM} + sigma_noise^2*eye(N_em(iM)));
            
            cons = [cons, (Hm{iM})'*SumJamming_approx_alpha*Hm{iM}>=10^(-20)];
            if (Method==0)
            cons = [cons, (Hm{iM})'*SumJamming_approx_beta*Hm{iM}>=10^(-20)];
            end

            cons = [cons, (Hm{iM})'*SumJamming_approx*Hm{iM}>=10^(-20)];
        else
            
            XI{2,iM} = sigma_noise^2*eye(N_em(iM));
        end
    else
    
    [Lambda1(iM), ~] = Get_Lambda({W1_current, V1_current, rho2p_current}, Hm{iM}, Gm2p{iM},0);
    [Lambda2(iM), ~] = Get_Lambda({W2_current, V2_current, rho1q_current}, Hm{iM}, Gm1q{iM},0);
    
    [am(iM), bm(iM), cm(iM), dm(iM)] = Get_abcd_m({W1_next, V1_next, rho2p_next}, {W2_next, V2_next, rho1q_next}, X1, X2, Hm{iM}, Gm2p{iM}, Gm1q{iM}, SCSI); % for f_mathcal and g_mathcal
    
    Lambda_approx(iM) = (2/alpha_current)*am(iM) + (2/beta_current)*bm(iM) - cm(iM)*alpha_next/(alpha_current^2)-dm(iM)*beta_next/(beta_current^2) + N_em(iM)*sigma_noise^2;
    
    if (~isempty(EH))
        disp('Add EH constraint hereeeee');
        EHm(iM) = delta*(Lambda_approx(iM) -  N_em(iM)*sigma_noise^2);
        if (~SumEH)
            cons = [cons, EHm(iM)>=Em(iM)];
        end
    end
    
    end
    
    zone = 1; % consider all users in zone 1
    
    for iK = 1:1:K
        
        if (WCS)
            
            gammahat_ED_1k(iM, iK) = (W1_current(:,iK))'*Hm{iM}*inv(XI{1,iM})*(Hm{iM})'*W1_current(:,iK);
            
%             Theta_ED_1k{iM, iK} = inv(XI{1,iM})*(Hm{iM})'*W1_current(:,iK)*W1_current(:,iK)'*Hm{iM}*inv(XI{1,iM});
            
%             C_hat_EDL_1k(iM, iK) = real(log(1+gammahat_ED_1k(iM,iK)) + (1+gammahat_ED_1k(iM,iK))^(-1)*( 2*(W1_current(:,iK))'*Hm{iM}*inv(XI{1,iM})*(Hm{iM})'*(W1_next(:,iK)- W1_current(:,iK)) ...
%                                                                                                         +trace(Theta_ED_1k{iM, iK}*XI{1,iM})-mu_m1k_next(iM, iK) ) );
            
            cons = [cons, cone([mu_m1k_next(iM, iK), 0.5*(Theta_ED_1k(iM,iK)-alpha_next)], 0.5*(Theta_ED_1k(iM,iK)+alpha_next)) ];
            
            cons = [cons, [(2*mu_m1k_next(iM, iK)*mu_m1k_current(iM, iK)-mu_m1k_current(iM, iK)^2), (W1_next(:,iK))'*Hm{iM}; (Hm{iM})'*W1_next(:,iK), ((Hm{iM})'*SumJamming_approx*Hm{iM}+ sigma_noise^2*eye(N_em(iM)) ) ] >= 0 ];
            
            cons = [cons, mu_m1k_next(iM, iK)>= 10^(-30)];
            
            cons = [cons, (2*mu_m1k_next(iM, iK)*mu_m1k_current(iM, iK)-mu_m1k_current(iM, iK)^2)>= 10^(-30)];
                
            C_hat_EDL_1k(iM, iK) = real(log(1+gammahat_ED_1k(iM,iK)) + (1+gammahat_ED_1k(iM,iK))^(-1)*( Theta_ED_1k(iM,iK) - gammahat_ED_1k(iM,iK) ) );
            
            cons = [cons, Theta_ED_1k(iM,iK)>=10^(-30)];
            
%             cons = [cons, Theta_ED_1k(iM,iK)<=10^(2)];
                                                                                                    
            cons = [cons, C_hat_EDL_1k(iM, iK) <= GammaDL_1k(iK)];

            cons = [cons, C_hat_EDL_1k(iM, iK)>=10^(-30)];     
            
        else
        
            % constraint (37a)

            if (SCSI)
                cons = [cons, cone([W1_next(:,iK)'*diag(sqrt(diag(Hm{iM}*Hm{iM}'))), 0.5*(beta_1k_DL(iK)-beta_bar_1k_DL(iM, iK)) ], 0.5*(beta_1k_DL(iK)+beta_bar_1k_DL(iM, iK)) ) ];
                cons = [cons, beta_bar_1k_DL(iM, iK)- mu_m1k_next(iM, iK) <= ((1-eps_m1k^(1/M))*N_em(iM)*sigma_noise^2)*alpha_next ];
                xi_beta_1k_DL(iK) = log(1+beta_1k_DL_current(iK)) + (1+beta_1k_DL_current(iK))^(-1)*(beta_1k_DL(iK)-beta_1k_DL_current(iK));
                cons = [cons, xi_beta_1k_DL(iK) <= GammaDL_1k(iK)];
            else

                gamma1k_EDL(iM,iK) = Get_gammaEDL(X_current, alpha_current, iK, Hm{iM}, sigma_noise, [Lambda1(iM) Lambda2(iM)], zone);

                cons = [cons, mu_p_m1k(iM, iK)>=10^(-30)];

                cons = [cons, cone([((Hm{iM})'*W1_next(:,iK))' 0.5*(mu_m1k_next(iM, iK) - mu_p_m1k(iM, iK))], 0.5*(mu_m1k_next(iM, iK) + mu_p_m1k(iM, iK))) ];

                C_tilde_EDL_1k(iM, iK) = real(log(1+gamma1k_EDL(iM,iK)) + (1+gamma1k_EDL(iM,iK))^(-1)*(mu_p_m1k(iM, iK) - gamma1k_EDL(iM,iK)));

                cons = [cons, C_tilde_EDL_1k(iM, iK) <= GammaDL_1k(iK)];

                cons = [cons, C_tilde_EDL_1k(iM, iK)>=10^(-30)];

            end
            
%         end

        
            % constraint (39)

            cons = [cons, mu_hat_m1k(iM, iK)>= 10^(-30)];


            cons = [cons, cone([ sqrt(alpha_p_current/mu_m1k_current(iM, iK))*mu_m1k_next(iM, iK), sqrt(mu_m1k_current(iM, iK)/alpha_p_current)*alpha_p, 0.5*(mu_hat_m1k(iM, iK)-1)], 0.5*(mu_hat_m1k(iM, iK)+1)) ];

%         if (WCS)
%              
%             F_mathcal_1k(iM, iK) = trace(Theta_ED_1k{iM, iK}*(Hm{iM}'*SumJamming_approx*Hm{iM} + sigma_noise^2*eye(N_em(iM)) ) );
%             
%         else
            
            if (SCSI)

                F_mathcal_1k(iM, iK) = Lambda_approx(iM) - (2/alpha_current)*real(W1_next(:,iK)'*diag(diag(Hm{iM}*Hm{iM}'))*W1_current(:,iK)) + real(W1_current(:,iK)'*diag(diag(Hm{iM}*Hm{iM}'))*W1_current(:,iK))*alpha_next/(alpha_current^2);

            else

                F_mathcal_1k(iM, iK) = Lambda_approx(iM) - (2/alpha_current)*real(W1_next(:,iK)'*Hm{iM}*Hm{iM}'*W1_current(:,iK)) + norm(Hm{iM}'*W1_current(:,iK))^2*alpha_next/(alpha_current^2);

            end
            
%         end
         
        cons = [cons, 0.5*mu_hat_m1k(iM, iK) <= F_mathcal_1k(iM, iK)];


        % constraint (44c)

        cons = [cons, mu_m1k_next(iM, iK)>= 10^(-30)];

%         cons = [cons, mu_m1k_next(iM, iK)<= 10^(4)];

        end
        
    end
    
    
    for iQ = 1:1:Q
        
        if (WCS)
            
            gammahat_EU_1q(iM, iQ) = (rho1q_current(iQ))^2*Gm1q{iM}(iQ,:)*inv(XI{2,iM})*(Gm1q{iM}(iQ,:))';
            
%             Theta_tilde_EU_1q{iM, iQ} = (rho1q_current(iQ))^2*inv(XI{1,iM})*(Gm1q{iM}(iQ,:))'*Gm1q{iM}(iQ,:)*inv(XI{1,iM});
            
%             C_hat_EUL_1q(iM, iQ) = real(log(1+gammahat_EU_1q(iM,iQ)) + (1+gammahat_EU_1q(iM,iQ))^(-1)*( 2*rho1q_current(iQ)*Gm1q{iM}(iQ,:)*inv(XI{1,iM})*(Gm1q{iM}(iQ,:))'*(rho1q_next(iQ)- rho1q_current(iQ)) ...
%                                                                                                         +trace(Theta_tilde_EU_1q{iM, iQ}*XI{1,iM})-mu_tilde_m1q_next(iM, iQ) ) );
                                                                                                    
%             cons = [cons, real(2*rho1q_current(iQ)*Gm1q{iM}(iQ,:)*inv(XI{1,iM})*(Gm1q{iM}(iQ,:))'*(rho1q_next(iQ)- rho1q_current(iQ)))>=10^(-20)];

            cons = [cons, cone([mu_tilde_m1q_next(iM, iQ), 0.5*(Theta_tilde_EU_1q(iM,iQ)-beta_next)], 0.5*(Theta_tilde_EU_1q(iM,iQ)+beta_next)) ];
            
            cons = [cons, [(2*mu_tilde_m1q_next(iM, iQ)*mu_tilde_m1q_current(iM, iQ)-mu_tilde_m1q_current(iM, iQ)^2), (rho1q_next(:,iQ)*Gm1q{iM}(iQ,:)); (rho1q_next(:,iQ)*(Gm1q{iM}(iQ,:))'), ((Hm{iM})'*SumJamming_approx*Hm{iM}+ sigma_noise^2*eye(N_em(iM)) ) ] >= 0 ];
            
            cons = [cons, mu_tilde_m1q_next(iM, iQ)>= 10^(-30)];
            
            cons = [cons, (2*mu_tilde_m1q_next(iM, iQ)*mu_tilde_m1q_current(iM, iQ)-mu_tilde_m1q_current(iM, iQ)^2)>= 10^(-30)];
            
            C_hat_EUL_1q(iM, iQ) = real(log(1+gammahat_EU_1q(iM,iQ)) + (1+gammahat_EU_1q(iM,iQ))^(-1)*( Theta_tilde_EU_1q(iM,iQ) - gammahat_EU_1q(iM,iQ) ) );
            
            cons = [cons, Theta_tilde_EU_1q(iM,iQ)>=10^(-30)];
            
%             cons = [cons, Theta_tilde_EU_1q(iM,iQ)<=10^(2)];
                                                                                                    
            cons = [cons, C_hat_EUL_1q(iM, iQ) <= GammaUL_1q(iQ)];

            cons = [cons, C_hat_EUL_1q(iM, iQ)>=10^(-30)];
            
        else
    
            % constraint (41b)

            if (SCSI)

                cons = [cons, cone([rho1q_next(iQ)*sqrt(Gm1q{iM}(iQ,:)*Gm1q{iM}(iQ,:)'), 0.5*(beta_1q_UL(iQ)-beta_bar_1q_UL(iM, iQ)) ], 0.5*(beta_1q_UL(iQ)+beta_bar_1q_UL(iM, iQ)) ) ];
                cons = [cons, beta_bar_1q_UL(iM, iQ)- mu_tilde_m1q_next(iM, iQ) <= ((1-eps_m1q^(1/M))*N_em(iM)*sigma_noise^2)*beta_next ];
                xi_beta_1q_UL(iQ) = log(1+beta_1q_UL_current(iQ)) + (1+beta_1q_UL_current(iQ))^(-1)*(beta_1q_UL(iQ)-beta_1q_UL_current(iQ));
                cons = [cons, xi_beta_1q_UL(iQ) <= GammaUL_1q(iQ)];

            else

                gamma1q_EUL(iM,iQ) = Get_gammaEUL(X_current, beta_current, iQ, Gm1q{iM}(iQ,:), sigma_noise, [Lambda1(iM) Lambda2(iM)], zone);

                cons = [cons, mu_tilde_p_m1q(iM, iQ)>=10^(-30)];

                cons = [cons, cone([(rho1q_next(iQ)*Gm1q{iM}(iQ,:)), 0.5*(mu_tilde_m1q_next(iM, iQ) - mu_tilde_p_m1q(iM, iQ))], 0.5*(mu_tilde_m1q_next(iM, iQ) + mu_tilde_p_m1q(iM, iQ))) ];

                C_tilde_EUL_1q(iM, iQ) = log(1+gamma1q_EUL(iM,iQ)) + (1+gamma1q_EUL(iM,iQ))^(-1)*(mu_tilde_p_m1q(iM, iQ) - gamma1q_EUL(iM,iQ));

                cons = [cons, C_tilde_EUL_1q(iM, iQ) <= GammaUL_1q(iQ)];

                cons = [cons, C_tilde_EUL_1q(iM, iQ)>=10^(-30)];

            end
        
%         end
        % constraint (42b)

        cons = [cons, (mu_hat_m1q(iM, iQ))>= 10^(-30)];


        cons = [cons, cone([ sqrt(beta_p_current/mu_tilde_m1q_current(iM, iQ))*mu_tilde_m1q_next(iM, iQ), sqrt(mu_tilde_m1q_current(iM, iQ)/beta_p_current)*beta_p, 0.5*(mu_hat_m1q(iM, iQ)-1)], 0.5*(mu_hat_m1q(iM, iQ)+1)) ];

%         if (WCS)
%             G_mathcal_1q(iM, iQ) = trace(Theta_tilde_EU_1q{iM, iQ}*(Hm{iM}'*SumJamming_approx*Hm{iM}+sigma_noise^2*eye(N_em(iM))) );
%         else
            if (SCSI)
                G_mathcal_1q(iM, iQ) = Lambda_approx(iM) - (2/beta_current)*(rho1q_current(iQ)*rho1q_next(iQ)*norm(Gm1q{iM}(iQ,:))^2) + ((rho1q_current(iQ)^2*norm(Gm1q{iM}(iQ,:))^2))*(beta_next/(beta_current^2));
            else
                G_mathcal_1q(iM, iQ) = Lambda_approx(iM) - (2/beta_current)*(rho1q_current(iQ)*rho1q_next(iQ)*norm(Gm1q{iM}(iQ,:))^2) + ((rho1q_current(iQ)^2*norm(Gm1q{iM}(iQ,:))^2))*(beta_next/(beta_current^2));
            end
        
%         end
        
        cons = [cons, 0.5*mu_hat_m1q(iM, iQ) <= G_mathcal_1q(iM, iQ)];
        

        
        % constraint (44c)
        
        cons = [cons, mu_tilde_m1q_next(iM, iQ)>= 10^(-30)];
        
        end
    end
    
    
    
    zone = 2; % consider all users in zone 1
    
    for iL = 1:1:L
    
        if (WCS)
            gammahat_ED_2l(iM, iL) = (W2_current(:,iL))'*Hm{iM}*inv(XI{2,iM})*(Hm{iM})'*W2_current(:,iL);
            
%             Theta_ED_2l{iM, iL} = inv(XI{2,iM})*(Hm{iM})'*W2_current(:,iL)*W2_current(:,iL)'*Hm{iM}*inv(XI{2,iM});
%             
%             C_hat_EDL_2l(iM, iL) = real(log(1+gammahat_ED_2l(iM,iL)) + (1+gammahat_ED_2l(iM,iL))^(-1)*( 2*(W2_current(:,iL))'*Hm{iM}*inv(XI{2,iM})*(Hm{iM})'*(W2_next(:,iL)- W2_current(:,iL)) ...
%                                                                                                         +trace(Theta_ED_2l{iM, iL}*XI{2,iM})-mu_m2l_next(iM, iL) ) );
            
            cons = [cons, cone([mu_m2l_next(iM, iL), 0.5*(Theta_ED_2l(iM,iL)-beta_next)], 0.5*(Theta_ED_2l(iM,iL)+beta_next)) ];
            
            cons = [cons, [(2*mu_m2l_next(iM, iL)*mu_m2l_current(iM, iL)-mu_m2l_current(iM, iL)^2), (W2_next(:,iL))'*Hm{iM}; (Hm{iM})'*W2_next(:,iL), ((Hm{iM})'*SumJamming_approx*Hm{iM}+ sigma_noise^2*eye(N_em(iM)) ) ] >= 0 ];
            
            cons = [cons, mu_m2l_next(iM, iL)>= 10^(-30)];
            
            cons = [cons, (2*mu_m2l_next(iM, iL)*mu_m2l_current(iM, iL)-mu_m2l_current(iM, iL)^2)>= 10^(-30)];
                
            C_hat_EDL_2l(iM, iL) = real(log(1+gammahat_ED_2l(iM,iL)) + (1+gammahat_ED_2l(iM,iL))^(-1)*( Theta_ED_2l(iM,iL) - gammahat_ED_2l(iM,iL) ) );
            
            cons = [cons, Theta_ED_2l(iM,iL)>=10^(-30)];
            
%             cons = [cons, Theta_ED_2l(iM,iL)<=10^(2)];
                                                                                                    
            cons = [cons, C_hat_EDL_2l(iM, iL) <= GammaDL_2l(iL)];

            cons = [cons, C_hat_EDL_2l(iM, iL)>=10^(-30)];
            
        else
            % constraint (41a)

            if (SCSI)
                cons = [cons, cone([W2_next(:,iL)'*diag(sqrt(diag(Hm{iM}*Hm{iM}'))), 0.5*(beta_2l_DL(iL)-beta_bar_2l_DL(iM, iL)) ], 0.5*(beta_2l_DL(iL)+beta_bar_2l_DL(iM, iL)) ) ];
                cons = [cons, beta_bar_2l_DL(iM, iL)- mu_m2l_next(iM, iL) <= ((1-eps_m2l^(1/M))*N_em(iM)*sigma_noise^2)*beta_next ];
                xi_beta_2l_DL(iL) = log(1+beta_2l_DL_current(iL)) + (1+beta_2l_DL_current(iL))^(-1)*(beta_2l_DL(iL)-beta_2l_DL_current(iL));
                cons = [cons, xi_beta_2l_DL(iL) <= GammaDL_2l(iL)];
            else

            gamma2l_EDL(iM,iL) = Get_gammaEDL(X_current, beta_current, iL, Hm{iM}, sigma_noise, [Lambda1(iM) Lambda2(iM)], zone);

            cons = [cons, mu_p_m2l(iM, iL)>=10^(-30)];

            cons = [cons, cone([((Hm{iM})'*W2_next(:,iL))' 0.5*(mu_m2l_next(iM, iL) - mu_p_m2l(iM, iL))], 0.5*(mu_m2l_next(iM, iL) + mu_p_m2l(iM, iL))) ];

            C_tilde_EDL_2l(iM, iL) = log(1+gamma2l_EDL(iM,iL)) + (1+gamma2l_EDL(iM,iL))^(-1)*(mu_p_m2l(iM, iL) - gamma2l_EDL(iM,iL));

            cons = [cons, C_tilde_EDL_2l(iM, iL) <= GammaDL_2l(iL)];

            cons = [cons, C_tilde_EDL_2l(iM, iL)>=10^(-30)];

            end
            
%         end
        
        % constraint (42a)

            cons = [cons, (mu_hat_m2l(iM, iL))>= 10^(-30)];


        cons = [cons, cone([ sqrt(beta_p_current/mu_m2l_current(iM, iL))*mu_m2l_next(iM, iL), sqrt(mu_m2l_current(iM, iL)/beta_p_current)*beta_p, 0.5*(mu_hat_m2l(iM, iL)-1)], 0.5*(mu_hat_m2l(iM, iL)+1)) ];
        
%         if (WCS)
%             
%             G_mathcal_2l(iM, iL) = trace(Theta_ED_2l{iM, iL}*(Hm{iM}'*SumJamming_approx*Hm{iM}+sigma_noise^2*eye(N_em(iM)) ) );
%             
%         else
            
            if (SCSI)

                G_mathcal_2l(iM, iL) = Lambda_approx(iM) - (2/beta_current)*real(W2_next(:,iL)'*diag(diag(Hm{iM}*Hm{iM}'))*W2_current(:,iL)) + real(W2_current(:,iL)'*diag(diag(Hm{iM}*Hm{iM}'))*W2_current(:,iL))*beta_next/(beta_current^2);

            else

                G_mathcal_2l(iM, iL) = Lambda_approx(iM) - (2/beta_current)*real(W2_next(:,iL)'*Hm{iM}*Hm{iM}'*W2_current(:,iL)) + norm(Hm{iM}'*W2_current(:,iL))^2*beta_next/(beta_current^2);
            end
            
%         end
        
        cons = [cons,  0.5*mu_hat_m2l(iM, iL) <= G_mathcal_2l(iM, iL)];


        
        % constraint (44c)
        
        cons = [cons, mu_m2l_next(iM, iL)>= 10^(-30)];
        
%         cons = [cons, mu_m2l_next(iM, iL)<= 10^(4)];
        end
        
    end
    
    
    
    for iP = 1:1:P
        
        if (WCS)
            gammahat_EU_2p(iM, iP) = (rho2p_current(iP))^2*Gm2p{iM}(iP,:)*inv(XI{1,iM})*(Gm2p{iM}(iP,:))';
            
%             Theta_tilde_EU_2p{iM, iP} = (rho2p_current(iP))^2*inv(XI{2,iM})*(Gm2p{iM}(iP,:))'*Gm2p{iM}(iP,:)*inv(XI{2,iM});
%             
%             C_hat_EUL_2p(iM, iP) = real(log(1+gammahat_EU_2p(iM,iP)) + (1+gammahat_EU_2p(iM,iP))^(-1)*( 2*rho2p_current(iP)*Gm2p{iM}(iP,:)*inv(XI{2,iM})*(Gm2p{iM}(iP,:))'*(rho2p_next(iP)- rho2p_current(iP)) ...
%                                                                                                         +trace(Theta_tilde_EU_2p{iM, iP}*XI{2,iM})-mu_tilde_m2p_next(iM, iP) ) );
            
                                                                                                    
            cons = [cons, cone([mu_tilde_m2p_next(iM, iP), 0.5*(Theta_tilde_EU_2p(iM,iP)-alpha_next)], 0.5*(Theta_tilde_EU_2p(iM,iP)+alpha_next)) ];
            
            cons = [cons, [(2*mu_tilde_m2p_next(iM, iP)*mu_tilde_m2p_current(iM, iP)-mu_tilde_m2p_current(iM, iP)^2), (rho2p_next(:,iP)*Gm2p{iM}(iP,:)); (rho2p_next(:,iP)*(Gm2p{iM}(iP,:))'), ((Hm{iM})'*SumJamming_approx*Hm{iM}+ sigma_noise^2*eye(N_em(iM)) ) ] >= 0 ];
            
            cons = [cons, mu_tilde_m2p_next(iM, iP)>= 10^(-30)];
            
            cons = [cons, (2*mu_tilde_m2p_next(iM, iP)*mu_tilde_m2p_current(iM, iP)-mu_tilde_m2p_current(iM, iP)^2)>= 10^(-30)];
            
            C_hat_EUL_2p(iM, iP) = real(log(1+gammahat_EU_2p(iM,iP)) + (1+gammahat_EU_2p(iM,iP))^(-1)*( Theta_tilde_EU_2p(iM,iP) - gammahat_EU_2p(iM,iP) ) );
            
            cons = [cons, Theta_tilde_EU_2p(iM,iP)>=10^(-30)];
            
%             cons = [cons, Theta_tilde_EU_2p(iM,iP)<=10^(1)];
            
            cons = [cons, C_hat_EUL_2p(iM, iP) <= GammaUL_2p(iP)];

            cons = [cons, C_hat_EUL_2p(iM, iP)>=10^(-30)];
        else
    
            % constraint (41c)
            if (SCSI)
                cons = [cons, cone([rho2p_next(iP)*sqrt(Gm2p{iM}(iP,:)*Gm2p{iM}(iP,:)'), 0.5*(beta_2p_UL(iP)-beta_bar_2p_UL(iM, iP)) ], 0.5*(beta_2p_UL(iP)+beta_bar_2p_UL(iM, iP)) ) ];
                cons = [cons, beta_bar_2p_UL(iM, iP)- mu_tilde_m2p_next(iM, iP) <= ((1-eps_m2p^(1/M))*N_em(iM)*sigma_noise^2)*alpha_next ];
                xi_beta_2p_UL(iP) = log(1+beta_2p_UL_current(iP)) + (1+beta_2p_UL_current(iP))^(-1)*(beta_2p_UL(iP)-beta_2p_UL_current(iP));
                cons = [cons, xi_beta_2p_UL(iP) <= GammaUL_2p(iP)];
            else
            gamma2p_EUL(iM,iP) = Get_gammaEUL(X_current, alpha_current, iP, Gm2p{iM}(iP,:), sigma_noise, [Lambda1(iM) Lambda2(iM)], zone);

            cons = [cons, mu_tilde_p_m2p(iM, iP)>=10^(-30)];

            cons = [cons, cone([(rho2p_next(iP)*Gm2p{iM}(iP,:)), 0.5*(mu_tilde_m2p_next(iM, iP) - mu_tilde_p_m2p(iM, iP))], 0.5*(mu_tilde_m2p_next(iM, iP) + mu_tilde_p_m2p(iM, iP))) ];

            C_tilde_EUL_2p(iM, iP) = log(1+gamma2p_EUL(iM,iP)) + (1+gamma2p_EUL(iM,iP))^(-1)*(mu_tilde_p_m2p(iM, iP) - gamma2p_EUL(iM,iP));

            cons = [cons, C_tilde_EUL_2p(iM, iP) <= GammaUL_2p(iP)];

            cons = [cons, C_tilde_EUL_2p(iM, iP)>=10^(-30)];

            end
%         end
        
        % constraint (42c)

            cons = [cons, (mu_hat_m2p(iM, iP))>= 10^(-30)];

        cons = [cons, cone([ sqrt(alpha_p_current/mu_tilde_m2p_current(iM, iP))*mu_tilde_m2p_next(iM, iP), sqrt(mu_tilde_m2p_current(iM, iP)/alpha_p_current)*alpha_p, 0.5*(mu_hat_m2p(iM, iP)-1)], 0.5*(mu_hat_m2p(iM, iP)+1)) ];
        
%         if (WCS)
%             F_mathcal_2p(iM, iP) = trace(Theta_tilde_EU_2p{iM, iP}*(Hm{iM}'*SumJamming_approx*Hm{iM}+sigma_noise^2*eye(N_em(iM))) );
%         else
            if (SCSI)
                F_mathcal_2p(iM, iP) = Lambda_approx(iM) - (2/alpha_current)*(rho2p_current(iP)*rho2p_next(iP)*norm(Gm2p{iM}(iP,:))^2) + ((rho2p_current(iP)^2*norm(Gm2p{iM}(iP,:))^2))*alpha_next/(alpha_current^2);
            else

                F_mathcal_2p(iM, iP) = Lambda_approx(iM) - (2/alpha_current)*(rho2p_current(iP)*rho2p_next(iP)*norm(Gm2p{iM}(iP,:))^2) + ((rho2p_current(iP)^2*norm(Gm2p{iM}(iP,:))^2))*alpha_next/(alpha_current^2);
            end
%         end
        
        cons = [cons, 0.5*mu_hat_m2p(iM, iP) <= F_mathcal_2p(iM, iP)];
        
        
        % constraint (44c)
        
        cons = [cons, mu_tilde_m2p_next(iM, iP)>= 10^(-30)];
        
        end
    end
    
    
    
end

if (SumEH)
    cons = [cons, sum(EHm)>=Em];
end

end


%% power constraints

% constraint (43a)

if (Method==3)
% disp('analyzing4 ..........................')    
%     Pbs_group1 = sdpvar(1, 1, 'full','real');
%     cons = [cons, cone([vec(W1_next); vec(V1_next); 0.5*(Pbs_group1-1)], 0.5*(Pbs_group1+1) ) ];
    if (IsDownlink)
    cons = [cons, cone([vec(W1_next); vec(V1_next)], sqrt(Pbs) )  ];
    end
%     cons = [cons, Pbs_group1<=Pbs];
% disp('analyzing5 ..........................')
else

    Pbs_group1 = sdpvar(1, 1, 'full','real');
    Pbs_group2 = sdpvar(1, 1, 'full','real');

    cons = [cons, cone([vec(W1_next); vec(V1_next); 0.5*(Pbs_group1-1)], 0.5*(Pbs_group1+1) ) ];
    cons = [cons, cone([vec(W2_next); vec(V2_next); 0.5*(Pbs_group2 - beta_next)], 0.5*(Pbs_group2 + beta_next)) ];

    Pbs_group1_appx = sdpvar(1, 1, 'full','real');

    for iK = 1:1:K
        if (iK==1)
            Pbs_group1_appx = real((W1_current(:,iK))'*W1_next(:,iK));
        else
            Pbs_group1_appx = Pbs_group1_appx + real((W1_current(:,iK))'*W1_next(:,iK));
        end
    end

    Pbs_group1_appx = Pbs_group1_appx + real(trace(V1_current'*V1_next));

    cons = [cons, (Pbs_group1+Pbs_group2+((norm(vec(W1_current)))^2 + (norm(vec(V1_current)))^2)*beta_next/(beta_current^2)-(2/beta_current)*Pbs_group1_appx) <= Pbs];

end

% constraint (43b)

if (Method~=3)

Pu_2p_appx = sdpvar(1,P,'full', 'real');

for iP = 1:1:P
    
    Pu_2p_appx(iP) = -(2/beta_current)*rho2p_current(iP)*rho2p_next(iP) + (rho2p_current(iP)^2)*(beta_next/(beta_current^2));
    
    cons = [cons, cone( [rho2p_next(iP) 0.5*(P_up_max(iP)-Pu_2p_appx(iP) - 1)], 0.5*(P_up_max(iP)-Pu_2p_appx(iP) + 1) ) ];
    
end

end

%% Optimization

% cons

if ((~isempty(findstr(preStatus{1},'Successfully'))) && preStatus{2}>0.1)
%     myops = sdpsettings('solver','mosek','verbose',0);
%     Status = preStatus;
    myops = sdpsettings('solver','sdpt3','verbose',0);
else
    myops = sdpsettings('solver','sdpt3','verbose',0);
end

% if (~isempty(EH) && ~isempty(findstr(preStatus{3},'NoEav')))
%     cons = [cons, eta>=0.5];
%     diagnotics = solvesdp(cons, -etaEH, myops)
% else
diagnotics = solvesdp(cons, -obj, myops)
% end
% diagnotics = solvemoment(cons, -obj, myops)

% if (~isempty(EH) && ~isempty(findstr(preStatus{3},'NoEav')))
%     OptimalValue = double(etaEH);
% else
    OptimalValue = double(eta);
% end




opt_X = {{double(W1_next), double(V1_next), double(rho2p_next)},{double(W2_next), double(V2_next), double(rho1q_next)}};
opt_alpha = double(alpha_next);
opt_beta = double(beta_next);
opt_mu = {double(mu_m1k_next), double(mu_m2l_next), double(mu_tilde_m1q_next), double(mu_tilde_m2p_next)};

C1k = double(C1k_DL)
C2l = double(C2l_DL)
C1q = double(C1q_UL)
C2p = double(C2p_UL)

Gamma_1k = double(GammaDL_1k)
Gamma_2l = double(GammaDL_2l)
Gamma_1q = double(GammaUL_1q)
Gamma_2p = double(GammaUL_2p)


% mu_m1k = double(mu_m1k_next)
% mu_m2l = double(mu_m2l_next)
% mu_tilde_m1q = double(mu_tilde_m1q_next)
% mu_tilde_m2p = double(mu_tilde_m2p_next)

alpha_p_next = double(alpha_p)
beta_p_next = double(beta_p)

% beta_dot_next = double(beta_dot)

if (M && isempty(findstr(preStatus{3},'NoEav')))
% Out_mu_p_m1k = double(mu_p_m1k)
% Out_mu_p_m2l = double(mu_p_m2l)
% Out_mu_p_m1q = double(mu_tilde_p_m1q)
% Out_mu_hat_m1k = double(mu_hat_m1k)
% Out_mu_hat_m2l = double(mu_hat_m2l)
% Out_mu_hat_m1q = double(mu_hat_m1q)
% Out_mu_hat_m2p = double(mu_hat_m2p)

if (WCS)
    Out_C_hat_EDL_1k = double(C_hat_EDL_1k)
    Out_C_hat_EDL_2l = double(C_hat_EDL_2l)
    Out_C_hat_EUL_1q = double(C_hat_EUL_1q)
    Out_C_hat_EUL_2p = double(C_hat_EUL_2p)
    
    Out_Theta_ED_1k = double(Theta_ED_1k)
    Out_Theta_tilde_EU_2p = double(Theta_tilde_EU_2p)
else
    Out_C_tilde_EDL_1k = double(C_tilde_EDL_1k)
    Out_C_tilde_EDL_2l = double(C_tilde_EDL_2l)
    Out_C_tilde_EUL_1q = double(C_tilde_EUL_1q)
    Out_C_tilde_EUL_2p = double(C_tilde_EUL_2p)
    if (~isempty(EH))
        
        Out_EHm = double(EHm)
        
    end
end

% Out_alpha_bar_DL = double(alpha_bar_DL)
% Out_alpha_bar_UL = double(alpha_bar_UL)
% Out_beta_bar_DL = double(beta_bar_DL)
% Out_beta_bar_UL = double(beta_bar_UL)
% 
% Out_alpha_hat = double(alpha_hat)
% Out_beta_hat = double(beta_hat)
% 
% Out_alpha_hat_p = double(alpha_hat_p)
% Out_beta_hat_p = double(beta_hat_p)

    if (SCSI)
        beta_SCSI_next = {[double(beta_1k_DL)], [double(beta_2l_DL)], [double(beta_1q_UL)], [double(beta_2p_UL)]};
    else
        beta_SCSI_next = {[],[],[],[]};
    end
    
else
    if (SCSI)
        beta_SCSI_next = {zeros(1,K),zeros(1,L),zeros(1,Q),zeros(1,P)};
    else
        beta_SCSI_next = {[],[],[],[]};
    end
end



Status = diagnotics.info;


end

