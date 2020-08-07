%%% Written by Van-Dinh Nguyen
%%% Authors: Van-Dinh Nguyen, Hieu V. Nguyen, Octavia A. Dobre, and Oh-Soon Shin
%%% Corresponding to the paper: "A New Design Paradigm for Secure Full-Duplex Multiuser Systems," 
%%%IEEE Journal on Selected Areas in Communications (SI in PLS for 5G Wireless Networks), vol. 36 , no. 7 , pp. 1480-1498, July 2018. 
%%%  




function [isBreak, Alg1_OptValue, Alg1_OptValueChain, OptSolution ] = ProposedAlg1( Pbs, P_up_max, AllChannels, sigma, MaxIteration, Fixed_Timegroup, Method, Ru, eps_SCSI, EH)
%PROPOSEDALG Summary of this function goes here
%   Detailed explanation goes here

if (nargin < (nargin('ProposedAlg1')-1))
    eps_SCSI = [];
    EH = [];
elseif (nargin < nargin('ProposedAlg1'))
    EH = [];
end

K = size(AllChannels{1, 1}, 2);
L = size(AllChannels{1, 2}, 2);
P = size(AllChannels{1, 3}, 2);
Q = size(AllChannels{1, 4}, 2);
M = size(AllChannels{1, 7}, 2);

N_t = size(AllChannels{1, 1}, 1);
N_r = size(AllChannels{1, 3}, 1);

if (Method==3)
    MaxIteration = 30;
end


strdisp = ['For SI = ' num2str(10*log10(sigma(1))) ' dB'];

if (~isempty(Fixed_Timegroup))
    strdisp = ['Fixed time -- ' strdisp];
    if (length(Fixed_Timegroup)==1)
        Fixed_Timegroup = [0.5 0.5];
    end
else
    strdisp = ['Opt. time -- ' strdisp];
end

global_OptValue = -1000;

    
        
        
        
    
        disp(['***** Getting a feasible point ..... ' strdisp]);

%         [ GetBreak, OptimalValue, OptimalValue_preStep, DownlinkRate_PerGroupPerUser, UplinkRate_PerGroupPerUser, RDown_current, RThDown_current, RUp_current, RTh_current, W_current, p_current, phi_current, time_current, w_tilde_current, p_bar_current, alpha_bar_current, alpha, beta_bar_current, beta] = GetInitialization(Pbs, P, Hu_tilde, Hd_tilde, G_SI_tilde, G_lk, W_current, p_current, phi_current, rho, Fixed_timegroup_assignment, Rate_Threshold);
    
%         [ GetBreak, OptimalValue, OptimalValue_preStep, UplinkRate_PerGroupPerUser, DownlinkRate_PerGroupPerUser, W_current, p_current, mu_current, vartheta_current] = GetInitialization(G, Pbs, P, Hu_tilde, Hd_tilde, G_SI_tilde, G_lk, error_lk, rho, Fixed_timegroup_assignment, Rate_Threshold, [biomega1 biomega2], eta);
        


        %% Loop for optimization
        
        GetBreak = 0;

        if (GetBreak)

            isBreak = 1;
            global_OptValue = 0;
            global_OptValueChain = [];
            global_OptValueChain_sum = 0;
            global_UplinkRate_PerGroupPerUser = 0;
            global_DownlinkRate_PerGroupPerUser = 0;
            
            

        else

            isBreak = 0;

            if (Method==3)
                
                disp(['---------------------------------- ' strdisp ' - Run iterative algorithm --------------------------------']);
                
                [DownRate, OptValueChain, X_current, alpha_current, beta_current, mu_current, time_next] = IterativeAlg(Pbs, P_up_max, AllChannels, sigma, MaxIteration, Fixed_Timegroup, [Method 1], N_t, K, 0, 0, 0, M, strdisp, Ru, eps_SCSI, EH);
                if (DownRate <=0)
                    UpRate = 0;
                    OptimalValue_current = 0;
                else
                [UpRate, ~, ~, ~, ~, ~, ~] = IterativeAlg(Pbs, P_up_max, AllChannels, sigma, MaxIteration, Fixed_Timegroup, [Method 0], N_t, 0, 0, 0, Q, M, strdisp, Ru, eps_SCSI);
                OptimalValue_current = min(DownRate,UpRate)/2;
                end
                disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> HD hereeeeeeeeeee');
                DownRate
                UpRate
                OptimalValue_current
                
                
            else            
                
                disp(['---------------------------------- ' strdisp ' - Run iterative algorithm --------------------------------']);
                
                [OptimalValue_current, OptValueChain, X_current, alpha_current, beta_current, mu_current, time_next] = IterativeAlg(Pbs, P_up_max, AllChannels, sigma, MaxIteration, Fixed_Timegroup, Method, N_t, K, L, P, Q, M, strdisp, Ru, eps_SCSI, EH);

            end
            
            disp('-------------------- Converging point ----------------------------');

            time_current = time_next

            Alg1_OptValue = 1/log(2)*OptimalValue_current                     
            Alg1_OptValueChain = OptValueChain/log(2)

            OptSolution = {X_current, alpha_current, beta_current, mu_current};

            disp('------------------------------------------------------------------');


        end



end

function [OptimalValue_current, OptValueChain, X_current, alpha_current, beta_current, mu_current, time_next] = IterativeAlg(Pbs, P_up_max, AllChannels, sigma, MaxIteration, Fixed_Timegroup, Method, N_t, K, L, P, Q, M, strdisp, Ru, eps_SCSI, EH)
%         eps_SCSI
        if (length(eps_SCSI)>1)
            SCSI = 1;
            WCS = 0;
        else
            if(length(eps_SCSI)==1)
                WCS=1;
            else
                WCS=0;
            end
            SCSI = 0;
        end
        
        if ((length(Method)==1) && (Method==0 && (10*log10(sigma(1))<-80)))
            Threshold = 1;
            MaxIteration = 50;
        else
            if (WCS)
                Threshold = 1.5;
            else
                if (~isempty(EH))
                    Threshold = 1;%-10^(10);
                else
                    Threshold = 2;
                end
            end
        end
        
        
            n = 0;
            n1= 0;
%             OptimalValue_current = OptimalValue;
            OptimalValue_current = 0;
            OptValueChain = [];
            
%             if (~isempty(EH))
%                 RunningMode = '';
%             else
                RunningMode = 'NoEav';
%             end
            preStatus = {'', 0, RunningMode};
            InitPoint = -1;
            
            W1_current = 20*(randn(N_t,K)+1i*randn(N_t,K)); %100
            W2_current = 20*(randn(N_t,L)+1i*randn(N_t,L));
            
            if (WCS)
                disp('>>>>>>>> WCS');
                V1_current = 10*(randn(N_t,N_t)+1i*randn(N_t,N_t));
                V2_current = 10*(randn(N_t,N_t)+1i*randn(N_t,N_t));

                rho1q_current = 20*rand(1, Q);
                rho2p_current = 20*rand(1, P);
            else
                
                if (~isempty(EH))
                    V1_current = 1*(randn(N_t,N_t)+1i*randn(N_t,N_t));
                    V2_current = 1*(randn(N_t,N_t)+1i*randn(N_t,N_t));
                else
                    if (Method==4)
                        V1_current = 0*(randn(N_t,N_t)+1i*randn(N_t,N_t));
                        V2_current = 0*(randn(N_t,N_t)+1i*randn(N_t,N_t));
                    else
                        V1_current = 1*(randn(N_t,N_t)+1i*randn(N_t,N_t));
                        V2_current = 1*(randn(N_t,N_t)+1i*randn(N_t,N_t));
                    end
                end

                rho1q_current = 3*rand(1, Q);
                rho2p_current = 3*rand(1, P);
            
            end
            
            X_current = {{W1_current, V1_current, rho2p_current},{W2_current, V2_current, rho1q_current}};
            
            if (~isempty(Fixed_Timegroup))
                alpha_current = 1/Fixed_Timegroup(1);
                beta_current = 1/Fixed_Timegroup(2);
            else
                alpha_current = 2;
                beta_current = 2;
            end
            
            if ((length(Method)==1) && (Method==0 && (10*log10(sigma(1))<-80)))
                [mu_m1k_current, mu_m2l_current, mu_tilde_m1q_current, mu_tilde_m2p_current] = Get_mu(X_current, alpha_current, beta_current, AllChannels, [10^(-7) sigma(2)], K, L, P, Q, M, WCS);
            else
                [mu_m1k_current, mu_m2l_current, mu_tilde_m1q_current, mu_tilde_m2p_current] = Get_mu(X_current, alpha_current, beta_current, AllChannels, sigma, K, L, P, Q, M, WCS);
            end
            
            mu_current = {mu_m1k_current, mu_m2l_current, mu_tilde_m1q_current, mu_tilde_m2p_current};

            alpha_p_current = 1/alpha_current;
            beta_p_current = 1/beta_current;
            
            if (SCSI)
                beta_SCSI_current = {[zeros(1,K)], [zeros(1,L)], [zeros(1,Q)], [zeros(1,P)]};
            else
                beta_SCSI_current = {[],[],[],[]};
            end
            
%             beta_dot_current = 1/(beta_current-1);
            OptimalValue = 0;

            while (n<MaxIteration)

                disp(['******************* ' strdisp '-- Iteration: ' num2str(n+1) '--' num2str(n1) ' *********************']);
                % Solve (44)

%                 [ OptimalValue, X_next, alpha_next, beta_next, mu_next,
%                 Status ] = Get_optSolutionPerIteration(G, Pbs, P_up_max, AllChannels, sigma, Fixed_Timegroup, X_current, alpha_current, beta_current, mu_current, preStatus); 
                if ((length(Method)==1) && (OptimalValue<1 && Method==0 && (10*log10(sigma(1))<-130)))
                    if (SCSI)
                        [ OptimalValue, X_next, alpha_next, beta_next, mu_next, alpha_p_next, beta_p_next, beta_SCSI_next,  Status ] = Get_optSolutionPerIteration2(Pbs, P_up_max, AllChannels, [10^(-7) sigma(2)], Fixed_Timegroup, X_current, alpha_current, beta_current, mu_current, alpha_p_current, beta_p_current, preStatus, Method, Ru, beta_SCSI_current, eps_SCSI, EH);
                    else
                        [ OptimalValue, X_next, alpha_next, beta_next, mu_next, alpha_p_next, beta_p_next, beta_SCSI_next, Status ] = Get_optSolutionPerIteration2(Pbs, P_up_max, AllChannels, [10^(-7) sigma(2)], Fixed_Timegroup, X_current, alpha_current, beta_current, mu_current, alpha_p_current, beta_p_current, preStatus, Method, Ru, beta_SCSI_current, eps_SCSI, EH);
                    end
                else
                    if (SCSI) 
                        [ OptimalValue, X_next, alpha_next, beta_next, mu_next, alpha_p_next, beta_p_next, beta_SCSI_next, Status ] = Get_optSolutionPerIteration2(Pbs, P_up_max, AllChannels, sigma, Fixed_Timegroup, X_current, alpha_current, beta_current, mu_current, alpha_p_current, beta_p_current, preStatus, Method, Ru, beta_SCSI_current, eps_SCSI, EH);
                    else
                        [ OptimalValue, X_next, alpha_next, beta_next, mu_next, alpha_p_next, beta_p_next, beta_SCSI_next, Status ] = Get_optSolutionPerIteration2(Pbs, P_up_max, AllChannels, sigma, Fixed_Timegroup, X_current, alpha_current, beta_current, mu_current, alpha_p_current, beta_p_current, preStatus, Method, Ru, beta_SCSI_current, eps_SCSI, EH);
                    end
                end
                
%                 if (OptimalValue==0 || (abs(OptimalValue)>10^(6) && n==5) || ((~isempty(findstr(RunningMode,'NoEav'))) && n==10))
%                     OptimalValue_current = 0;
%                     break;
%                 end

                if (OptimalValue==0)
                    OptimalValue_current = 0;
                    break;
                end

                % Update
                disp('hereeeeeeeeeeeeeeeee');
                time_next = 1./[alpha_next beta_next]
                OptimalValue_current
                OptimalValue
                
                if (isFeasible(time_next, OptimalValue_current, OptimalValue, InitPoint))


                    X_current = X_next;
                                     
                    alpha_current = alpha_next
                    
                    beta_current = beta_next
                    
                    W1W1 = (1/alpha_current)*trace(X_current{1,1}{1,1}'*X_current{1,1}{1,1})
                    W2W2 = (1/beta_current)*trace(X_current{1,2}{1,1}'*X_current{1,2}{1,1})
                    
                    V1V1 = (1/alpha_current)*trace(X_current{1,1}{1,2}'*X_current{1,1}{1,2})
                    V2V2 = (1/beta_current)*trace(X_current{1,2}{1,2}'*X_current{1,2}{1,2})
                    
                    SumPbs = W1W1 + V1V1 + W2W2 + V2V2
                    
                    rho2p_sqr = (1/alpha_current)*X_current{1,1}{1,3}.^2
                    rho1q_sqr = (1/beta_current)*X_current{1,2}{1,3}.^2
                    
                    if (isempty(findstr(preStatus{3},'NoEav')))
                    
                        mu_current = mu_next;
                        
                    else
                        
                        [mu_m1k_current, mu_m2l_current, mu_tilde_m1q_current, mu_tilde_m2p_current] = Get_mu(X_current, alpha_current, beta_current, AllChannels, sigma, K, L, P, Q, M, WCS);
            
                        mu_current = {mu_m1k_current, mu_m2l_current, mu_tilde_m1q_current, mu_tilde_m2p_current};
                        
                    end
                    
                    alpha_p_current = alpha_p_next;
                    beta_p_current = beta_p_next;
%                     beta_dot_current = beta_dot_next;

                    if (SCSI)
                        beta_SCSI_current = beta_SCSI_next;
                    end

                    if ((~isempty(findstr(Status,'Successfully'))) && OptimalValue>Threshold && (~isempty(findstr(RunningMode,'NoEav'))))
                        if (length(Method)>1)
                            X_current{1,1}{1,2} = 50*(randn(N_t,N_t)+1i*randn(N_t,N_t));
                        end
                        RunningMode = '';                     
                    end
                    
                    if ((Method==2) && OptimalValue>Threshold && (~isempty(findstr(RunningMode,'NoEav'))))
                        if (length(Method)>1)
                            X_current{1,1}{1,2} = 50*(randn(N_t,N_t)+1i*randn(N_t,N_t));
                        end
                        RunningMode = '';                     
                    end
                    
                    if (isempty(findstr(RunningMode,'NoEav')))
                        InitPoint = InitPoint+1;
                    end

                    preStatus = {Status, OptimalValue, RunningMode};


                    % Check convergence

        %             if (abs(OptimalValue-OptimalValue_current)<10^-2)
                    if (checkConvergence(OptValueChain, OptimalValue, n) && (OptimalValue>0))
                        break;
                    else
                        
                        if (OptimalValue>0 && OptimalValue>=OptimalValue_current)
                            OptValueChain = [OptValueChain OptimalValue];
                        else
                            if (OptimalValue>2 && OptimalValue<OptimalValue_current)
                                break
                            end
                            OptValueChain = [];
                        end
                        OptimalValue_current = OptimalValue
%                         OptValueChain_sum = [OptValueChain_sum (sum(sum(UplinkRate_PerGroupPerUser)) + sum(sum(DownlinkRate_PerGroupPerUser)))];
                    end

                else

                    disp('--> keeping the latest feasible point');
                    OptValueChain = [OptValueChain OptimalValue_current];
                    OptimalValue = OptimalValue_current;
%                     OptValueChain_sum = [OptValueChain_sum (sum(sum(UplinkRate_PerGroupPerUser)) + sum(sum(DownlinkRate_PerGroupPerUser)))];
                    break;

                end
                
                if (OptimalValue<0)
                    n1 = n1+1;
                    if (n1>50)
                        break
                    end
                else
                    n = n + 1;
                end

            end
            
    if ((OptimalValue<0) || (~isempty(findstr(RunningMode,'NoEav'))))
        OptimalValue_current = 0;
    end
    
%     if (OptimalValue<Threshold)
%         OptimalValue_current = 0;
%     end
end

function [check] = isFeasible(time_next, OptimalValue_current, OptimalValue, InitPoint) %time, , UplinkRate, DownlinkRate, Threshold

% check if the current point is feasible and the next point is infeasible

check = 1;

% if ((sum(( (UplinkRate-Threshold) >= -10^(-1) )) < length(UplinkRate)) || (sum(( (UplinkRate) < 100 )) < length(UplinkRate)))
%     check = 0;
%     disp('Uplink --> Not pass Threshold check');
% end
% 
% if ((sum(( (DownlinkRate-Threshold) >= -10^(-1) )) < length(DownlinkRate)) || (sum(( (DownlinkRate) < 100 )) < length(DownlinkRate)))
%     check = 0;
%     disp('Downlink --> Not pass Threshold check');
% end
% check
% check_time = (round(sum(time),3)<=1.000)
check_time =1;

% if (sum(time_next>1)>0 || sum(time_next<0)>0)
%     check_time = 0;
%     disp('Time_next --> Not pass range check');
% end

% check1 = [(sum(check_time) == 0]; % check1 = 0 means that current point is feasible


% check_time_next = (round(sum(time_next),3)<=1.000)
check_time_next = 1;

% check2 = [(sum(sum(check_alpha_next)) + sum(sum(check_beta_next)) + check_time_next) == 0]; % check2 = 1 means that next point is feasible

check3 = 1;

% if ((OptimalValue_current>OptimalValue && OptimalValue_current>0 && InitPoint>0) || OptimalValue==0)
if ((OptimalValue_current>2 && OptimalValue<0 && InitPoint>0) || OptimalValue==0)
    check3 = 0;
end

% if ((OptimalValue_current>0 && OptimalValue<0) || (round(OptimalValue_current,20)==0) || (OptimalValue_current>1000))
%     check3 = 0;
%     disp('OptimalValue --> Not pass range check');
%     OptimalValue_current
%     OptimalValue
% end

check = (check_time && check_time_next && check3 && check);

end


function [check] = checkConvergence(OptValueChain, OptimalValue, n)

check = 0;

% if (length(OptValueChain)>30 && (n>30 || (n>10 && OptimalValue>1) ) )
%     if (abs(OptimalValue-OptValueChain(length(OptValueChain))) < 10^-3)
%         check = 1;
%     end
% %     if (abs(OptimalValue-OptValueChain(length(OptValueChain)))/abs(OptimalValue) < 10^-2)
% %         check = 1;
% %     end
%     if (abs(OptimalValue-OptValueChain(length(OptValueChain)-5)) < 0.01)
%         check = 1;
%     end
% end

end

