addpath('../UtilityFunctions');

%% - flags

make_weights = 1;
simulate_network = 1;

pattern_completion_curve_sim = 1;
pattern_completion_curve_plot = 1;

plot_results = 1;    

%% - network params
dt = 1;
tau = 10;

NE = 400;
NI = 400;
N = NE+NI;

conn_spar_EE = 1;

N_itr = 30;
eta = 1/5;

m = 0;
mEE = m;
mEI = m/2;
mIE = m/2;
mII = m/2;

po_exc = linspace(0,pi,NE);
po_inh = linspace(0,pi,NI);
po_all = [po_exc, po_inh];

pert_size = .1;

l0 = 0.8;

ts_dur = [10,20,50,100];
N_rep = 10;
t_betw = 50; 

Ns_pert =[10,20,50,100,200];

% - sample plot
k1_s = 2;
k2_s = 3;

%%

r_all = {};
w_all = {};
evm = zeros(2,N_itr);
evecm = zeros(2, N);
evecmall = zeros(2,N_itr+1,N);
ccm = {};
ccm0 = {};

ks = [1, 4];
for gi = [2, 1]
k = ks(gi)

J0 = 1/NE*2;
JEE = J0 /conn_spar_EE;
JEI = J0 *k;
JIE = -J0 *k;
JII = -J0 *k;

%% - weight matrix
if make_weights

    wEE = zeros(NE,NE);
    wEI = zeros(NE,NI);
    for i = 1:NE
        wEE(i,:) = (1 + mEE * cos(2*(po_exc(i) - po_exc))) * JEE;
        wEI(i,:) = (1 + mEI * cos(2*(po_exc(i) - po_inh))) * JEI;
    end
    wIE = zeros(NI,NE);
    wII = zeros(NI,NI);
    for i = 1:NI
        wIE(i,:) = (1 + mIE * cos(2*(po_inh(i) - po_exc))) * JIE;
        wII(i,:) = (1 + mII * cos(2*(po_inh(i) - po_inh))) * JII;
    end

    wEE = wEE .* binornd(1,conn_spar_EE,[NE,NE]);
    
    w = [wEE, wEI
         wIE, wII];

    w = w .*(1+1*(rand(N,N)-.5));

    w(1:NE,:) = rectify(w(1:NE,:));
    w(NE+1:end,:) = -rectify(-w(NE+1:end,:));

    w(eye(N)==1) = 0;
end
    
%% - simulate
if simulate_network
   
    N_pert = Ns_pert(k1_s);
    t_dur = ts_dur(k2_s);
    t_betw = t_dur;
    
    wp = w;
        
    [V,D] = eig(w);
    ev_before = diag(D);
    [evh,b] = max(real(diag(D)));
    evecmall(gi,1,:) = V(:,b);
    
    % - perturb.
    pert_ids = 1:N_pert;
    
    if k == 1; pert_size = .1;
    elseif k == 4; pert_size = .1;
    end
    
    w_all{gi}{1} = w;
    evm(gi,1) = max(real(ev_before));
    for kk = 1:N_itr
        
        [r,v, dw] = pert_sim_short(N, NE, wp, dt, tau, pert_ids, pert_size, t_dur, N_rep, t_betw);

        r_all{gi}{kk} = r;

        zz = dw(1:NE,1:N);
        zz(isnan(zz))=0;

        wp_test = wp;
        wp_test(1:NE,1:NE) = wp_test(1:NE,1:NE) + eta * zz(1:NE,1:NE);
        
        [V,D] = eig(wp_test);
        evh = max(real(diag(D)));
        if evh < l0
            wp = wp_test;
        end

        w_all{gi}{kk+1} = wp;
        [V,D] = eig(wp);
        [evh,b] = max(real(diag(D)));
        evm(gi, kk+1) = evh;
        evecmall(gi,kk+1,:) = V(:,b);
    end
    
    [V,D] = eig(wp);
    [a,b] = max(real(diag(D)));
    evecm(gi,:) = V(:,b);
    
    % - pattern completion
    
    pert_ids = 1:N_pert/2;
        
    N_rep = 1;
    t_dur = 300;
    t_betw = 0;
    
    if k == 1; pert_size = 1;
    elseif k == 4; pert_size = .1;
    end
    
    % - just one
    pert_ids = 1:N_pert/2;
    
    % before induction/learning
    [r_pc0{gi},~, ~] = pert_sim_short(N, NE, w, dt, tau, pert_ids, pert_size, t_dur, N_rep, t_betw);

    % after induction/learning
    [r_pc{gi},~, ~] = pert_sim_short(N, NE, wp, dt, tau, pert_ids, pert_size, t_dur, N_rep, t_betw);
            
    % - for different fractions
    if pattern_completion_curve_sim
        ccm{gi} = zeros(1,N_pert-1);
        ccm0{gi} = zeros(1,N_pert-1);
        for i = 1:N_pert-1
            pert_ids = 1:i;
            [r,~, ~] = pert_sim_short(N, NE, wp, dt, tau, pert_ids, pert_size, t_dur, N_rep, t_betw);
            
            cc = corr(r(1:NE,:)');
            
            r0 = nanmean(r(:,300:300+t_dur),2);
            rp = nanmean(r(:,50:300),2);
            dr = (rp - r0); 
            
            ccm{gi}(i) = nanmean(dr(i+1:N_pert)) ./ nanmean(dr(1:i)); 
            ccm0{gi}(i) = nanmean(dr(N_pert+1:NE)) ./ nanmean(dr(1:i)); 
        end
    end
end

end

%% - 

if plot_results
    
    for gi = 1:2
        k = ks(gi)
        
    % - sample weights
    ws_id = [0, N_itr/2, N_itr]+1;
    ttls = {'Initial', 'Middle', 'Final'};
    
    figure('Position',[100,100,500,140])
    for i = 1:3
    subplot(1,3,i); 
    if k == 1
        title(ttls{i}, 'FontWeight', 'normal');
    end
    hold on
    
    imagesc(w_all{gi}{ws_id(i)}(1:N_pert,1:N_pert));
    colorbar()
    axis image
    
    xticks([1,N_pert]);
    yticks([1,N_pert]);
    
    if i ~= 1
        xticklabels([]); yticklabels([]);
    end
    if i == 3
        xlabel('post #');
        ylabel('pre #');
    end
    
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out', 'ydir', 'normal')
    end
    print(['k' num2str(k) '_sampleW.png'], '-dpng', '-r300');
    
    % - pattern completion
    figure('Position',[100,100,270,225])
    
    hold on
    plot(r_pc{gi}(1+N_pert:NE,:)', '-', 'color', [.7,.7,.7]);
    plot(r_pc{gi}(1:N_pert/2,:)', '-', 'color', 'r');
    plot(r_pc{gi}(1+N_pert/2:N_pert,:)', '-', 'color', [1,0.5,0]);
    
    xlabel('Time')
    ylabel('Activity')
    
    if k == 4
        yticks([0, .1, .2, .3])
        ylim([0, .3]);
    else
        ylim([0,8])
    end
    
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out', 'ydir', 'normal')
    
    print(['k' num2str(k) '_patternCompl.png'], '-dpng', '-r300');
    
    % - before induction:
    if k == 1
        figure('Position',[100,100,250,225])
    elseif k == 4
        figure('Position',[100,100,270,225])
    end
    
    hold on
    plot(r_pc0{gi}(1+N_pert:NE,:)', '-', 'color', [.7,.7,.7]);
    plot(r_pc0{gi}(1:N_pert/2,:)', '-', 'color', 'r');
    plot(r_pc0{gi}(1+N_pert/2:N_pert,:)', '-', 'color', [1,0.5,0]);
    
    xlabel('Time')
    ylabel('Activity')
    
    if k == 4
        yticks([0, .1, .2, .3])
        ylim([0, .3]);
    elseif k == 1
        ylim([0,8]);
    end
    
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out', 'ydir', 'normal')
    
    print(['k' num2str(k) '_patternCompl_beforeInduction.png'], '-dpng', '-r300');
    
    % - pattern completion curve
    if pattern_completion_curve_plot
    
    figure('Position',[100,100,225,225])
    hold on
    
    h1 = plot((1:N_pert-1)/N_pert*100, ccm{gi}*100, '-o', 'linewidth',2, 'color', [1,0.5,0]);
    h2 = plot((1:N_pert-1)/N_pert*100, ccm0{gi}*100, '-o', 'linewidth',2, 'color', [.7,.7,.7]);
    
    if k == 4
        h3 = plot((1:N_pert-1)/N_pert*100, ccm{1}*100, '--', 'linewidth',2, 'color', [1,0.5,0]);
        h4 = plot((1:N_pert-1)/N_pert*100, ccm0{1}*100, '--', 'linewidth',2, 'color', [.7,.7,.7]);
    end
    
    if k == 4
        legend([h1,h2], {'within', 'outside'}, 'location', 'best')
        legend boxoff
    end
    
    xlabel('Partial activ. (%)');
    ylabel('Fraction resp. (%)');
    
    yticks([0, 25, 50, 75, 100])
    
    ylim([-10,80])
    
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out')
    
    print(['k' num2str(k) '_patternCompl_curve.png'], '-dpng', '-r300');
    
    end
    
    end
    
end
