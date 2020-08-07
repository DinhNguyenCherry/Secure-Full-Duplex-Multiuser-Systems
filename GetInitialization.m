function [ X_next, alpha, beta ] = GetInitialization( Pbs, P_up_max, X_current, alpha_current, beta_current, K, L, P, Q, N_t )
%GETINITIALIZATION1 Summary of this function goes here
%   Detailed explanation goes here
Pbs 
P_up_max
%% extract current solution

X1 = X_current{1,1};
X2 = X_current{1,2};

W1_current = X1{1,1};
V1_current = X1{1,2};
rho2p_current = X1{1,3};

W2_current = X2{1,1};
V2_current = X2{1,2};
rho1q_current = X2{1,3};

%% initialize variables

W1_next = sdpvar(N_t, K, 'full','complex');
V1_next = sdpvar(N_t, N_t, 'full','complex');%zeros(N_t,N_t);%
rho2p_next = sdpvar(1, P, 'full','real');

W2_next = sdpvar(N_t, L, 'full','complex');
V2_next = sdpvar(N_t, N_t, 'full','complex');%zeros(N_t,N_t);%
rho1q_next = sdpvar(1, Q, 'full','real');

alpha_next = sdpvar(1, 1, 'full','real');
beta_next = sdpvar(1, 1, 'full','real');




%% Setting
obj = sum(rho2p_next);
cons = [];

%% Time constraints

alpha_p = sdpvar(1, 1, 'full','real');
beta_p = sdpvar(1, 1, 'full','real');

cons = [cons, alpha_next>=1];
cons = [cons, beta_next>=1];
cons = [cons, alpha_next<=100];
cons = [cons, beta_next<=100];

cons = [cons, cone([1, 0.5*(alpha_next - alpha_p)], 0.5*(alpha_next + alpha_p))];
cons = [cons, cone([1, 0.5*(beta_next - beta_p)], 0.5*(beta_next + beta_p))];
cons = [cons, alpha_p>=10^(-2)];
cons = [cons, beta_p>=10^(-2)];

cons = [cons, (alpha_p + beta_p) <= 1 ];

%% power constraints

for iQ = 1:1:Q
    cons = [cons, rho1q_next(iQ) >= 10^(-4)];
end

for iP = 1:1:P
    cons = [cons, rho2p_next(iP) >= 10^(-4)];
end

for iQ = 1:1:Q    
    cons = [cons, cone([rho1q_next(iQ), 0.5*(P_up_max(P+iQ)-beta_next)], 0.5*(P_up_max(P+iQ)+beta_next)) ];   
end

% constraint (43a)

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


% constraint (43b)

Pu_2p_appx = sdpvar(1,P,'full', 'real');

for iP = 1:1:P
    
    Pu_2p_appx(iP) = -(2/beta_current)*rho2p_current(iP)*rho2p_next(iP) + (rho2p_current(iP)^2)*(beta_next/(beta_current^2));
    
    cons = [cons, cone( [rho2p_next(iP) 0.5*(P_up_max(iP)-Pu_2p_appx(iP) - 1)], 0.5*(P_up_max(iP)-Pu_2p_appx(iP) + 1) ) ];
    
end


myops = sdpsettings('solver','sdpt3','verbose',0);
diagnotics = solvesdp(cons, -obj, myops)


X_next = {{double(W1_next), double(V1_next), double(rho2p_next)},{double(W2_next), double(V2_next), double(rho1q_next)}};
alpha = 1/(1-1/double(beta_next));
beta = double(beta_next);

end

