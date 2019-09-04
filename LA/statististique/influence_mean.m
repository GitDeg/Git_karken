%% Two-sample t test
close all
iterations = 50; 

for mu = 1 : 10
    
    fig = figure (mu);
    for suj = 1 : 12
        
        n_suj = 1 : 12;
        
        Y1 = ALL_mean_1_5_manu(mu+1).grp2(1:12, 1:end-1);
        %Y2 = ALL_mean_1_5_manu(mu+1).grp1(:, 1:end-1);
        Y2(n_suj(n_suj~=suj),:) = ALL_mean_1_5_manu(mu+1).grp2(n_suj(n_suj~=suj), 1:end-1 );
        Y2(suj,:) = mean(ALL_mean_1_5_manu(mu+1).grp2(n_suj(n_suj~=suj), 1:end-1 ));
        
        A = [ones(size(Y1,1),1); 2*ones(size(Y2,1),1)];
        
        t = spm1d.stats.nonparam.anova1([Y1;Y2], A); %,'equal_var',false
        ti = t.inference(0.05, 'iterations', iterations);
        
        
        subplot(6,2,suj)
        hold on
        spm1d.plot.plot_meanSD(Y1, 'color' ,'r');
        spm1d.plot.plot_meanSD(Y2, 'color', 'g');
        
        ti.plot()
        
        for i = 1 : size(ti.clusters,2)
            SP = round(ti.clusters{1, i}.endpoints);
            line([SP(1),SP(1)],ylim,'Color',[1 0 0 0.1])
            line([SP(2),SP(2)],ylim,'Color',[1 0 0 0.1])
        end
        
    end
        
end