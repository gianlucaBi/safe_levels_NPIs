% This function implements the functional relationship in the dynamics.
% Namely, given:
%       \dot x = f(t,x),
% it evaluates f(t,x)


function dz = dynamics(t,z,m,params)
    aux_params      = num2cell(params);
    [beta, gamma, eps, delta, sigma, k_ih, k_hd, rho, nu, eta, k_id, y, vax_uptake, N, Ncomp, h_cap]  = deal(aux_params{:});
    
    gain        = 1;
    
    x           = z(1:Ncomp*N);
    u           = z(Ncomp*N+1:(Ncomp+1)*N);
    zResh       = reshape(x,Ncomp,N); 
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
    
    
    
    % Controller 1 Dynamics
    du  = zeros(N,1);
    for kk=1:N
        Lu     = 2*(u(kk)-1); 
        pp      = u(kk) - gain*Lu;
        pp      = min(pp,(.98*h_cap(kk)-m(2,kk))/m(1,kk));
        if h(kk) > h_cap(kk)           % when constraint is violated set u=0
            pp      = min(pp,0);
        end
        pp          = max(pp,0);    % project onto interval [0,1]
        pp          = min(pp,1);    
        du(kk)      = pp - u(kk);
    end

                

    
    dz      = [reshape([ds'; de'; di'; dh'; dr'; dv'; dd],Ncomp*N,1); du];
end