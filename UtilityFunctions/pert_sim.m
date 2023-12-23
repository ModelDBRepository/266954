function [r,v, dw, avg_ind, avg_ind_all, tot_ind, tot_ind_all, avg_ind_w] = ...
    pert_sim(N, NE, w, dt, tau, pert_ids, pert_size, t_dur, N_rep, t_betw, learning_rule)
   
if nargin < 11
    learning_rule = 'covariance';
end
    
    t_trans = 50;
    t0 = 300;

    N_pert = length(pert_ids);

    t_dur_tot = (t_dur+t_betw)*N_rep;
    
    t_end = t0 + t_dur_tot +t0;
    T = 0:dt:t_end;

    t_ids = [];
    for i = 1:N_rep
        t1 = t0 + (t_dur+t_betw)*(i-1);
        t2 = t1 + t_dur;
        t_ids = [t_ids, find((T>t1).*(T<t2))];
    end
    
    I_genPert = 1 + .1*(rand(N,length(T)));

    dI_pert = zeros(N,length(T));
    for i = 1:N_pert
        dI_pert(pert_ids(i), t_ids) = pert_size;
    end

    [r, v] = simulate_dynamics(I_genPert + dI_pert, N, T, w, dt, tau);

    t_avg = (t0:t_end-t0)/dt;
    t_base = (t_trans:t0)/dt;
    
    % - pre r x post v(/r)
    if strcmp(learning_rule, 'covariance') % - (pre - m) x (post - m)
    z1 = r(:,t_avg)-nanmean(r(:,t_base),2);
    z2 = v(:,t_avg)-nanmean(v(:,t_base),2);
  
    elseif strcmp(learning_rule, 'pre') % - (pre - m) x (post)
    z1 = r(:,t_avg)-nanmean(r(:,t_base),2);
    z2 = v(:,t_avg);
    
    elseif strcmp(learning_rule, 'post') % - (pre) x (post - m)
    z1 = r(:,t_avg);
    z2 = v(:,t_avg)-nanmean(v(:,t_base),2);
    end
    
    dw = (z1) * (z2)' / length(t_avg);
    
    % - correlation-based    
    dw(eye(N)==1)=nan;
    
    avg_ind = nanmean(nanmean(dw(pert_ids,pert_ids)));
    avg_ind_all = nanmean(nanmean(dw(1:NE,1:NE)));
    
    tot_ind = nanmean(nansum(dw(pert_ids,pert_ids)));
    tot_ind_all = nanmean(nansum(dw(1:NE,1:NE)));
    
    % -- weight-based
    L = (eye(N) - w);
    
    s = zeros(N,1);
    s(pert_ids) = pert_size;
    
    r_w = L\s;
    cc_w = r_w * r_w';
    cc_w(eye(N)==1)=nan;
    avg_ind_w = nanmean(nanmean(cc_w(pert_ids, pert_ids))) *t_dur/(t_dur+t_betw);
    
end

