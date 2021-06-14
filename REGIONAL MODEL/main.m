% Author: Gianluca Bianchin
% Date: June 14
% This code implements a primal-dual controller for the SEIHRS network 
% model with hospitalization + loss of immunity with N zones
%
% In summary, it seeks to chose u such that: (i) a capacity constraint on 
% the infections is satisfied at the peak, and (ii) a given cost is
% minimized at the peak. The relationship between u and the value of the
% peak is evaluated through linear regression



clear all
close all
clc


% Initializations
% global beta eps delta iC A iStar N sigma gammaH gammaR rhoR rhoD nu eta Ncomp phi
% Initial state
load interactionMatrix.mat
N           = size(A,1);                % Number of zones
Npop        = sum(regionPopSize);       % population size
gamma       = 1/9;
beta        = .58;                      
k_ih        = 0.02029464;                % Fraction of hospitalized 0.02029464
k_id        = 0.002162289;               % death rate outside hospitals 0.002162289
k_hd        = 0.099204;                 % death rate in hospitals 0.099204;
rho         = 1/8;                      % 1/hospitalization period
eps         = 1/4.2;                    % 1/incubation period
nu          = 0.9;                      % vaccine efficacy
sigma       = 1/(2*365);                % 1/natural immunity duration
eta         = 1/(2*365);                % 1/vaccination immunity duration 
vax_uptake  = .7;
delta       = 0.02965/365;              % death/birth rate
Ncomp       = 7;                        % number of compartments
state_hCap  = 8;                        % state hospitalization capacity, per 100K
h_cap       = state_hCap/100000*ones(11,1);
y           = 20000/Npop;               % .* regionPopSize'/sum(regionPopSize);                   % vaccination rate
weights     = ones(N,1);                % regional controller weights
params      = [beta gamma eps delta  sigma k_ih k_hd rho nu eta k_id vax_uptake  N Ncomp];




% Initial conditions of the combined population
load initCond.mat   % loads x0_March1 (initial state March 1)
x0          = x0_March1;
u0          = (1-.79)*ones(N,1);                                    % Initial state controller; (1-0.78) .22
x0          = reshape(x0,Ncomp*N,1); 
z0          = [x0; u0];



% Runge-Kutta parameters
tS      = 0;                % Initial time
tF      = 365;              % Final time
M       = 1*tF;              % Number of steps
h       = (tF-tS)/M;        % Step Size
tSpan   = tS:h:tF;          % Define time vector




%% Solution via Runge-Kutta
xRK         = zeros((Ncomp+1)*N,M);
xRK(:,1)    = z0;
for j=1:M     
        %% Estimate input to peak relationship
        x0Sim       = xRK(1:Ncomp*N,j);
        u_0         = xRK(Ncomp*N+1:(Ncomp+1)*N,j);   
        m           = zeros(2,N);
        for kk=1:N         % Cycle over the regions
            u_kVec              = u_0(kk) + [-.1:.05:.1];           % Perturb u(kk)
            u_kVec(u_kVec<=0)   = [];   u_kVec(u_kVec>1)   = [];    % Remove entries if negative or larger than 1
            u_kVec              = u_kVec - u_0(kk);
            
            nU                  = length(u_kVec);
            Ztop                = zeros(kk-1,nU);
            Zbottom             = zeros(N-kk,nU);
            uVec                = u_0 + [Ztop; u_kVec; Zbottom];
            
            X   = zeros(nU,2); Y = zeros(nU,1);
            for ii=1:nU     % Cycle over the inputs
                uSim            = uVec(:,ii);
                [tODE,zODE]     = ode23s(@(t,z) simulation_dynamics(t,z,uSim,params,A,h_cap,y),[0 15],x0Sim);       % Solve ODE 

                X(ii,:)         = [uVec(kk,ii) 1];
                iBar            = zODE(end,Ncomp*kk-3);
                Y(ii,1)         = iBar;
            end
            mAux        = pinv(X)*Y;        % iMax = m(1)*u+ m(2)
            m(:,kk)     = mAux;             % m(1,kk) contains m1 for zone kk

        end
        
        
        
        
        %% SEIR + Controller update step
        k1 = h*feval('dynamics', tSpan(j), xRK(:,j), m,params,A, h_cap,y,weights);
        k2 = h*feval('dynamics', tSpan(j)+h/2, xRK(:,j)+k1/2, m,params,A,h_cap,y,weights);
        k3 = h*feval('dynamics', tSpan(j)+h/2, xRK(:,j)+k2/2, m,params,A,h_cap,y,weights);
        k4 = h*feval('dynamics', tSpan(j)+h/2, xRK(:,j)+k3, m,params,A,h_cap,y,weights);
         
        xRK(:,j+1)   = xRK(:,j) + 1/6*(k1+2*k2+2*k3+k4);
end





%% Plots
S       = xRK(1:Ncomp:N*Ncomp,:);
E       = xRK(2:Ncomp:N*Ncomp,:);
I       = xRK(3:Ncomp:N*Ncomp,:);
H       = xRK(4:Ncomp:N*Ncomp,:);
R       = xRK(5:Ncomp:N*Ncomp,:);
V       = xRK(6:Ncomp:N*Ncomp,:);
D       = xRK(7:Ncomp:N*Ncomp,:);
U       = xRK(N*Ncomp+1:N*Ncomp+N,:);

myDates     = datetime(2021,3,1):days(1):datetime(2021,3,1)+tF;   % time span in date/time format
    

save workspace_LPHA_control

% Plot full-Colorado state
weighs = regionPopSize/sum(regionPopSize);
S_state = weighs*S;
E_state = weighs*E;
I_state = weighs*I;
H_state = weighs*H;
R_state = weighs*R;
V_state = weighs*V;
D_state = weighs*D;
U_state = weighs*U;
figure; 
subplot(8,1,1); hold on; plot(myDates,S_state,'linewidth',2); ylabel('s');
title('COMBINED COLORADO STATE')
subplot(8,1,2); hold on; plot(myDates,E_state,'linewidth',2); ylabel('e');
subplot(8,1,3); hold on; plot(myDates,I_state,'linewidth',2); ylabel('i');
subplot(8,1,4); hold on; plot(myDates,H_state,'linewidth',2); ylabel('h');
hold on; plot(myDates,ones(1,length(myDates)).*sum(h_cap.*regionPopSize')/Npop,'linewidth',2); 
subplot(8,1,5); hold on; plot(myDates,R_state,'linewidth',2); ylabel('r');
subplot(8,1,6); hold on; plot(myDates,V_state,'linewidth',2); ylabel('v');
subplot(8,1,7); hold on; plot(myDates,D_state,'linewidth',2); ylabel('v');
subplot(8,1,8); hold on; plot(myDates,U_state,'linewidth',2); ylabel('u');




%% Plot a single zone
plotZone    = 4;        % Desired zone you would like to plot
figure('name',strcat('States of ZONE = ',num2str(plotZone)))
subplot(8,1,1); hold on; plot(myDates,S(plotZone,:),'linewidth',2); ylabel('s');
title(strcat('States of ZONE ',num2str(plotZone)))
subplot(8,1,2); hold on; plot(myDates,E(plotZone,:),'linewidth',2); ylabel('e');
subplot(8,1,3); hold on; plot(myDates,I(plotZone,:),'linewidth',2); ylabel('i');
subplot(8,1,4); hold on; plot(myDates,H(plotZone,:),'linewidth',2); ylabel('h');
hold on; plot(myDates,ones(1,length(myDates)).*h_cap(plotZone),'linewidth',2);
subplot(8,1,5); hold on; plot(myDates,R(plotZone,:),'linewidth',2); ylabel('r');
subplot(8,1,6); hold on; plot(myDates,V(plotZone,:),'linewidth',2); ylabel('v');
subplot(8,1,7); hold on; plot(myDates,D(plotZone,:),'linewidth',2); ylabel('d');
subplot(8,1,8); hold on; plot(myDates,U(plotZone,:),'linewidth',2); ylabel('u');








%% Functions
function [dz,y] = simulation_dynamics(t,z,u,params,A,h_cap,y)
    aux_params      = num2cell(params);
    [beta, gamma, eps, delta, sigma, k_ih, k_hd, rho, nu, eta, k_id, vax_uptake, N, Ncomp]  = deal(aux_params{:});

    zResh       = reshape(z,Ncomp,N);
    s           = zResh(1,:)';
    e           = zResh(2,:)';
    i           = zResh(3,:)';
    h           = zResh(4,:)';
    r           = zResh(5,:)';
    v           = zResh(6,:)';
    
    theta       = 1-.23;                    % probability of vaccinating an individual in s compartment. Note: pVr+pVs=1
    ds          = - beta*diag(s)*A*diag(u)*i - theta*nu*y - delta*s + delta*ones(N,1) + sigma*r + eta*v;
    de          = - (delta+eps)*e + beta*diag(s)*A*diag(u)*i;
    di          = - (delta+gamma)*i + eps*e;
    dh          = - rho*h + k_ih*gamma*i;
    dr          = - (sigma+delta)*r - (1-theta)*nu*y + (1-k_ih-k_id)*gamma*i + (1-k_hd)*rho*h;
    dv          = - eta*v - delta*v + nu*y;
    dd          = k_id*gamma*i + k_hd*rho*h;
   
    dz          = reshape([ds'; de'; di'; dh'; dr'; dv'; dd'],Ncomp*N,1);
end


