function [r,v, dw] = pert_sim_short(N, NE, w, dt, tau, pert_ids, pert_size, t_dur, N_rep, t_betw, t0,I_inh)
    
if nargin < 12
I_inh=0;
end
if nargin < 11
    t0 = 300;
end

    t_trans = 50;
    
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
    
    I_genPert = 1 + .1/2*rand(N,length(T));
    I_genPert(NE+1:end,:) = I_genPert(NE+1:end,:) + I_inh;

    dI_pert = zeros(N,length(T));
    for i = 1:N_pert
        dI_pert(pert_ids(i), t_ids) = pert_size;
    end

    [r, v] = simulate_dynamics(I_genPert + dI_pert, N, T, w, dt, tau);

    t_avg = (t0:t_end-t0)/dt;
    t_base = (t_trans:t0)/dt;
    
    % - pre r x post v(/r)
    z1 = r(:,t_avg)-nanmean(r(:,t_base),2);
    z2 = v(:,t_avg)-nanmean(v(:,t_base),2);
    
    dw = (z1) * (z2') / length(t_avg);
    
    dw(eye(N)==1)=nan;
    
end