% Author: Gianluca Bianchin
% Date: March 12, 2021
% This script computes the number of days required to achieve u=1 at
% different levels of hospitalization capacity

clear all; close all; clc;

load('initCond')      
d0_Mar1             = d0;                                       % Inital number of deaths on Mar 1


h_cap_vec_100K      = 4:2:20;                               % Hospitalization capacity in ppl/100K
h_cap_vec           = h_cap_vec_100K/100000*5840795;        % Hospitalization capacity in %
y_vec               = [15000 20000 25000];                  % Vaccination rate
vax_uptk_vec        = .4:.05:.9;                            % Vaccination uptake
n_h                 = length(h_cap_vec);
n_y                 = length(y_vec);
n_u                 = length(vax_uptk_vec);
tF                  = 4*365;                                % Simulation horizon


%% Loops
values      = [];
for ii=1:n_h
    for jj=1:n_y
        for kk=1:n_u
            h_cap           = h_cap_vec(ii);
            y               = y_vec(jj);
            vax_uptk        = vax_uptk_vec(kk);
            workspaceName   = strcat('h_cap=',strrep(num2str(h_cap),'.','_'),'_y=',strrep(num2str(y),'.','_'),'v_uptk=',strrep(num2str(vax_uptk),'.','_'));

            [daysToFree, d_Aug1, tSpan]      = main_CO_1_zone(h_cap,y,vax_uptk,workspaceName,tF);

            aux             = [h_cap; y; vax_uptk; daysToFree; d_Aug1];
            values          = [values aux];         % EACH COLUMN OF VALUES:    [hosp_capacity; vax_rate; vax_uptk; days_to_free; deaths_Aug1]
        end
    end
end
save myWorkspace



%% Plots
fontSize = 18;
my_vax_uptk     = vax_uptk_vec(3);
for jj=1:n_y
    my_vax_rate    = y_vec(jj);
    ind = find(values(3,:)==my_vax_uptk  & values(2,:)==my_vax_rate);       % find indices with given vax_uptake and given vax_rate
    
    cap_vec             = values(1,ind);
    days_free_vec       = values(4,ind);
    deaths_vec          = values(5,ind);
    figure(1); hold on;  plot(cap_vec,days_free_vec,'linewidth',2,'DisplayName',strcat(num2str(y_vec(jj)),'vax/day'));
    figure(2); hold on;  plot(cap_vec,deaths_vec,'linewidth',2,'DisplayName',strcat(num2str(y_vec(jj)),'vax/day'));
end



%% Outputs
figure(1); ylabel('Days to $u=1$','FontUnits','points','interpreter','latex','FontWeight','normal','FontSize',fontSize,'FontName','Times'); legend('show');
figure(2); ylabel('Deaths August 1','FontUnits','points','interpreter','latex','FontWeight','normal','FontSize',fontSize,'FontName','Times'); legend('show');
 
