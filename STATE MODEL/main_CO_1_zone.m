function [daysToNormal, d_Aug1, tSpan] = main_CO_1_zone(hCap_p,y_p,v_uptake_p,workspaceName, tF)
% Date: March 10


Npop        = 5840795;                  % population size
gamma       = 1/9;
lambda      = 1.39;
beta        = .44*lambda;                      
k_ih        = 0.0143762;                % Fraction of hospitalized 0.02029464
k_id        = 0.0007142;               % death rate outside hospitals 0.002162289
k_hd        = 0.0421;                 % death rate in hospitals 0.099204;
rho         = 1/7.489;                  % 1/hospitalization period
eps         = 1/4.2;                    % 1/incubation period
delta       = 0.02965/365;              % death/birth rate
y           = y_p/Npop;                 % vaccination rate 
nu          = 0.81;                      % vaccine efficacy
sigma       = 1/(2*365);                % 1/natural immunity duration
eta         = 1/(2*365);                % 1/vaccination immunity duration
vax_uptake  = v_uptake_p;               % vaccination uptake
N           = 1;                        % Number of zones
Ncomp       = 7;                        % number of compartments
h_cap       = hCap_p/Npop;             % Hospitalization capacity (1077/Npop) 5250 1342
params      = [beta gamma eps delta  sigma k_ih k_hd rho nu eta k_id y vax_uptake  N Ncomp h_cap];

load('initCond')                        % loads (s0 e0 i0 r0 h0 v0) generated through script in 3_SEIHRS/1 ZONE/main_CO_1zone
x0          = kron(ones(1,N),[s0; e0; i0; h0; r0; v0; d0]);
u0          = (1-.79)*ones(N,1);                                    % Initial state controller; (1-0.78) .22
x0          = reshape(x0,Ncomp*N,1); 
z0          = [x0; u0];




% Runge-Kutta parameters
tS          = 0;                % Initial time
M           = 1*tF;             % Number of steps
h           = (tF-tS)/M;        % Step Size
tSpan       = tS:h:tF;          % Define time vector




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
                [tODE,zODE]     = ode45(@(t,z) simulation_dynamics(t,z,uSim,params),[0 15],x0Sim);       % Solve ODE 

                X(ii,:)         = [uVec(kk,ii) 1];
                iBar            = zODE(end,Ncomp*kk-3);
                Y(ii,1)         = iBar;
            end
            mAux        = pinv(X)*Y;        % iMax = m(1)*u+ m(2)
            m(:,kk)     = mAux;             % m(1,kk) contains m1 for zone kk
            mPl(:,j)    = m;             

        end
        
        
        
        
        %% SEIR + Controller update step
        k1 = h*feval('dynamics', tSpan(j), xRK(:,j), m,params);
        k2 = h*feval('dynamics', tSpan(j)+h/2, xRK(:,j)+k1/2, m,params);
        k3 = h*feval('dynamics', tSpan(j)+h/2, xRK(:,j)+k2/2, m,params);
        k4 = h*feval('dynamics', tSpan(j)+h/2, xRK(:,j)+k3, m,params);
         
        xRK(:,j+1)   = xRK(:,j) + 1/6*(k1+2*k2+2*k3+k4);
end



%% Plots
s           = xRK(Ncomp-6,:);
e           = xRK(Ncomp-5,:);
i           = xRK(Ncomp-4,:);
h           = xRK(Ncomp-3,:);
r           = xRK(Ncomp-2,:);
v           = xRK(Ncomp-1,:);
d           = xRK(Ncomp,:);
u           = xRK(Ncomp*N+1,:);
% 
% myDates     = datetime(2021,3,1):days(1):datetime(2021,3,1)+tF;   % time span in date/time format
% figure; 
% subplot(4,1,1); hold on
% plot(myDates,s,'linewidth',2);  plot(myDates,e,'linewidth',2);  plot(myDates,i,'linewidth',2); 
% plot(myDates,h,'linewidth',2);  plot(myDates,r,'linewidth',2); plot(myDates,v,'linewidth',2);
% legend({'$S$','$E$','$I$','$H$','$R$','$V$'},'interpreter','latex')
% 
% subplot(4,1,2)
% plot(myDates,u,'linewidth',2)
% legend({'$u(t)$'},'interpreter','latex')
% 
% subplot(4,1,3); hold on
% plot(myDates,Npop*h,'linewidth',2); 
% plot(myDates, Npop*h_cap*ones(1,length(tSpan)));
% legend({'$H$','$i_\textup{cap}$'},'interpreter','latex')
% 
% subplot(4,1,4); hold on
% plot(myDates,Npop*d,'linewidth',2); 
% legend({'$d$'},'interpreter','latex')
% % 



d_Aug1       = Npop*d(153);               % August 1 = March 1 + 153 days
% Determine days to u=1
u           = xRK(Ncomp*N+1,:);
indx        = find(abs(u-1)<1e-3);
if isempty(indx)
    daysToNormal      = inf;
else
    daysToNormal      = tSpan(indx(1));
end

save(workspaceName)

end



%% Functions
function [dz,y] = simulation_dynamics(t,z,u,params)
    aux_params      = num2cell(params);
    [beta, gamma, eps, delta, sigma, k_ih, k_hd, rho, nu, eta, k_id, y, vax_uptake, N, Ncomp, h_cap]  = deal(aux_params{:});

    zResh       = reshape(z,Ncomp,N);
    s           = zResh(1,:)';
    e           = zResh(2,:)';
    i           = zResh(3,:)';
    h           = zResh(4,:)';
    r           = zResh(5,:)';
    v           = zResh(6,:)';
    
    if v > vax_uptake       % vaccination uptake
        y=0;
    end
    
    
    theta       = 1-.23;                    % probability of vaccinating an individual in s compartment. Note: pVr+pVs=1
    ds          = - beta*diag(s)*diag(u)*i - theta*nu*y - delta*s + delta + sigma*r + eta*v;
    de          = - (delta+eps)*e + beta*diag(s)*diag(u)*i;
    di          = - (delta+gamma)*i + eps*e;
    dh          = - rho*h + k_ih*gamma*i;
    dr          = - (sigma+delta)*r - (1-theta)*nu*y + (1-k_ih-k_id)*gamma*i + (1-k_hd)*rho*h;
    dv          = - eta*v - delta*v + nu*y;
    dd          = k_id*gamma*i + k_hd*rho*h;
   
    dz          = reshape([ds'; de'; di'; dh'; dr'; dv'; dd],Ncomp*N,1);
end


