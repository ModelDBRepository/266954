addpath('../UtilityFunctions');

%% - flags

make_weights = 1;
simulate_network = 1;

plot_induction_curves = 1;
plot_activity_samples = 1;
plot_induction_specificity = 1;

%% - network params

dt = 1;
tau = 10;

NE = 500;
NI = 500;
N = NE+NI;

weight_spec = 1;
conn_spar_EE = 1;

if weight_spec == 1
    m = 0.5;
else
    m = 0;
end

mEE = m;
mEI = m;
mIE = m;
mII = m;

po_exc = linspace(0,pi,NE);
po_inh = linspace(0,pi,NI);
po_all = [po_exc, po_inh];

pert_size = .1;

ts_dur = [10,20,50,100];
N_rep = 10;

Ns_pert =[10,50,100,150,200];

% - sample plot of activity
k1_s = 2;
k2_s = 3;

%% - 

for k = [1, 4]
k

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

    w = w .*(1+1/2*(rand(N,N)-.5));

    w(1:NE,:) = rectify(w(1:NE,:));
    w(NE+1:end,:) = -rectify(-w(NE+1:end,:));

    w(eye(N)==1) = 0;
end

%% - simulate
if simulate_network
    
    r_all = [];
    v_all = [];
    dw_all = [];
    pert_ids_all = [];
    
    avg_ind = zeros(length(Ns_pert),length(ts_dur)); 
    avg_ind_all = zeros(size(avg_ind)); 
    tot_ind = zeros(length(Ns_pert),length(ts_dur)); 
    tot_ind_all = zeros(size(tot_ind)); 
    avg_ind_w = zeros(length(Ns_pert),length(ts_dur)); 
    
    for k1 = 1:length(Ns_pert)
        N_pert = Ns_pert(k1)

        for k2 = 1:length(ts_dur)
            t_dur = ts_dur(k2);
            t_betw = t_dur;
            
            % - perturbed neurons
            pert_ids = 1:NE;
            iz = 200;
            pert_ids = pert_ids(iz:iz+N_pert-1);
            pert_ids_all{k1,k2} = pert_ids;
            
            [r,v, dw, avg_ind(k1,k2), avg_ind_all(k1,k2), tot_ind(k1,k2), tot_ind_all(k1,k2), avg_ind_w(k1,k2)] = ...
             pert_sim(N, NE, w, dt, tau, pert_ids, pert_size, t_dur, N_rep, t_betw);
     
            r_all{k1,k2} = r;
            v_all{k1,k2} = v;
            dw_all{k1,k2} = dw;

        end
    end
end

%% - plot 
cls = {'r', 'k'};
    
if plot_induction_curves
for ki = 1:2
    if ki == 1
        ki_pre = 'avg';
    elseif ki == 2
        ki_pre = 'tot';
    end

    figure('Position', [100,100,225,225])
    subplot(111); hold on;
    if ki == 1
        zz = avg_ind;
    elseif ki == 2
        zz = tot_ind;
    end

    plot(Ns_pert, zz./nanmax(zz(:)), '-o', 'LineWidth',2)

    if k== 1
    legend(string(ts_dur), 'Location', 'northwest', 'FontSize', 12)
    legend boxoff
    end

    ylim([0,1]);
    xlabel('Ensemble size')

    if ki == 1
    ylabel('Avg. Pot. (norm.)')
    elseif ki == 2
    ylabel('Ensemble Pot. (norm.)')    
    end

    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out')

    print(['k' num2str(k) '_' ki_pre 'Induction.png'], '-dpng', '-r300');

end
end
    
% - 

if plot_activity_samples
    r = r_all{k1_s,k2_s};
    v = v_all{k1_s,k2_s};
    dw = dw_all{k1_s,k2_s};
    pert_ids = pert_ids_all{k1_s,k2_s};
    N_pert = length(pert_ids);

    figure('Position', [100,100,400,400])
    zz = dw(1:N,1:N);
    zz(isnan(zz))=0;
    zm = nanmax(abs(zz(:)));
    imagesc(zz, [-zm,zm]/2)
    xticks([1,400,800])
    yticks([1,400,800])
    colorbar('Location', 'northoutside', 'FontSize',15);
    colormap('redblue')
    xlabel('Post-syn. Neuron #'); 
    ylabel('Pre-syn. Neuron #');
    axis image
    set(gca, 'LineWidth', 1, 'FontSize', 20, 'Box', 'off', 'TickDir', 'out', 'YDir', 'normal');
    print(['k' num2str(k) '_dw.png'], '-dpng', '-r300');

    % - 
    figure('Position', [100,100,600,300])
    subplot(2,1,1); hold on
    imagesc(r(1+NE:end,:));
    hcb=colorbar();
    title(hcb,'Activity');
    xlim([100,1500]);
    ylim([1,NI]);
    yticks([1, NI])
    yticklabels([NE+1, N])
    xticklabels([]);
    ylabel('Inh. Neuron #')
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out');

    subplot(2,1,2); hold on
    imagesc(r(1:NE,:));
    colorbar()
    xlim([100,1500]);
    ylim([0,NE]);
    yticks([1, NE])
    yticklabels([1, NE])
    ylabel('Exc. Neuron #')
    xlabel('Time (ms)')
    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out');

    print(['k' num2str(k) '_sampleAct.png'], '-dpng', '-r300');

end    

% - 

if plot_induction_specificity

    N_itr = 100;
    tot_ind_rand = zeros(size(avg_ind));
    tot_ind_rand_s = zeros(size(tot_ind_rand));
    for k1 = 1:length(Ns_pert)
        dw = dw_all{k1,k2_s};
        pert_ids = pert_ids_all{k1,k2_s};
        N_pert = Ns_pert(k1);

        zz_r = zeros(1,N_itr);
        for ir = 1:N_itr
            pert_id_rand = randperm(NE);
            pert_id_rand = pert_id_rand(1:N_pert);

            zz_r(ir) = nanmean(nansum(dw(pert_ids,pert_id_rand)));
        end
        tot_ind_rand(k1,k2_s) = nanmean(zz_r(:));
        tot_ind_rand_s(k1,k2_s) = nanstd(zz_r(:));
    end

    figure('Position',[100,100,300,500])

    subplot(211); 
    hold on

    plot(Ns_pert, squeeze(tot_ind(:,k2_s)), '-o', 'color', [1,.5,0], 'LineWidth',2)
    errorbar(Ns_pert, squeeze(tot_ind_rand(:,k2_s)), squeeze(tot_ind_rand_s(:,k2_s)), ...
        '-', 'Color', [.7,.7,.7], 'LineWidth',2);

    xlabel('Ensemble size')
    ylabel('Ensemble Pot.')

    legend({'within-assemb.', 'out-of-assemb.'}, 'location','best')
    legend boxoff

    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out');

    subplot(212)
    xw1 = squeeze(tot_ind(:,k2_s));
    xa1 = squeeze(tot_ind_rand(:,k2_s));
    ind_ind1 = (xw1-xa1)./(xw1+xa1);

    bar(Ns_pert, ind_ind1, 'k');

    xlabel('Ensemble size')
    ylabel('Ensemble Spec.')
    
    ylim([0,1.1]);

    set(gca, 'LineWidth', 1, 'FontSize', 15, 'Box', 'off', 'TickDir', 'out');

    print(['k' num2str(k) '_inductSpec.png'], '-dpng', '-r300');
end

end