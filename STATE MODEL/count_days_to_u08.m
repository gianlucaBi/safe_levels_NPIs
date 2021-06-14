clear  all
close all
clc


load ('myWorkspace','n_h','n_y','n_u','h_cap_vec','y_vec','vax_uptk_vec')
values = [];


for ii=1:n_h
    for jj=1:n_y
        for kk=1:n_u
            h_cap           = h_cap_vec(ii);
            y               = y_vec(jj);
            vax_uptk        = vax_uptk_vec(kk);
            
            workspaceName   = strcat('h_cap=',strrep(num2str(h_cap),'.','_'),'_y=',strrep(num2str(y),'.','_'),'v_uptk=',strrep(num2str(vax_uptk),'.','_'));
            load(workspaceName,'Npop','d','xRK','N','tSpan','Ncomp')
            
            
            d_Aug1          = Npop*d(153);               % August 1 = March 1 + 153 days
            u               = xRK(Ncomp*N+1,:);
            indx            = find(abs(u-.8)<1e-2);
            indx(indx<50)   =[];                            % resolve numerical issue of peak 
            if isempty(indx)
                daysToNormal      = inf;
            else
                daysToNormal      = tSpan(indx(1));
            end
            
            figure; plot(u);

            aux             = [h_cap; y; vax_uptk; daysToNormal; d_Aug1];
            values          = [values aux];         % EACH COLUMN OF VALUES:    [hosp_capacity; vax_rate; vax_uptk; days_to_free; deaths_Aug1]
        end
    end
end


save myWorkspace_u08 values