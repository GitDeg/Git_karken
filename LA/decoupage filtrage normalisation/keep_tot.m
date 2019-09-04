

for mu = 1 : 10
    All_mean_grp(1).mu(mu).data = ALL_mean_1_5_manu(mu+1).grp1;
    All_mean_grp(1).mu(mu).name = ALL_mean_1_5_manu(mu+1).Names;
    
    All_mean_grp(2).mu(mu).data = ALL_mean_1_5_manu(mu+1).grp2;
    All_mean_grp(2).mu(mu).name = ALL_mean_1_5_manu(mu+1).Names;
    All_mean_grp(3).mu(mu).data = ALL_mean_1_5_manu(mu+1).grp3;
    All_mean_grp(3).mu(mu).name = ALL_mean_1_5_manu(mu+1).Names;
    All_mean_grp(4).mu(mu).data = ALL_mean_1_5_manu(mu+1).grp4;
    All_mean_grp(4).mu(mu).name = ALL_mean_1_5_manu(mu+1).Names;
    All_mean_grp(5).mu(mu).data = ALL_mean_1_5_manu(mu+1).grp5;
    All_mean_grp(5).mu(mu).name = ALL_mean_1_5_manu(mu+1).Names;
    All_mean_grp(6).mu(mu).data = ALL_mean_1_5_manu(mu+1).grp6;
    All_mean_grp(6).mu(mu).name = ALL_mean_1_5_manu(mu+1).Names;
    All_mean_grp(7).mu(mu).data = ALL_mean_1_5_manu(mu+1).grp7;
    All_mean_grp(7).mu(mu).name = ALL_mean_1_5_manu(mu+1).Names;
    All_mean_grp(8).mu(mu).data = ALL_mean_1_5_manu(mu+1).grp8;
    All_mean_grp(8).mu(mu).name = ALL_mean_1_5_manu(mu+1).Names;
end

%% ttt1 mini 35% de reste 

p = double(perte > 0.35);

for grp = 1 : 8
    for mu = 1 : 10
        
        data = All_mean_grp(grp).mu(mu).data(:,1:end-1);
        sd = std(data);
        
        figure(grp)
        hold on
        subplot(5,2,mu)
        
        plot(data','k')
        hold on 
        plot(mean(data), 'g', 'linewidth' , 2)
        plot(mean(data)+3*sd, 'r')
        
        p_i = []; 
        p_i = p(mu, grp, :);
        p_suj = 1 : 12;
        
        data = All_mean_grp(grp).mu(mu).data(p_suj(p_i(:)==1),1:end-1);
        sd = std(data);
        
        figure(grp)
        hold on
        subplot(5,2,mu)
        
        plot(data','w')
        plot(All_mean_grp(grp).mu(mu).data(p_suj(p_i(:)==0),1:end-1)','k')
        hold on 
        plot(mean(data), 'g', 'linewidth' , 2)
        plot(mean(data)+3*sd, 'r')
        legend(num2str(length(p_suj(p_i(:)==1))))
        
        All_mean_grp(grp).mu(mu).ttt1 = All_mean_grp(grp).mu(mu).data;
        All_mean_grp(grp).mu(mu).ttt1(p_suj(p_i(:)==0), 1:end-1) = NaN;
        
    end
end

%% ttt2 influence sur la variance 

for grp = 1 : 8
    for mu = 1 : 10
        
        All_mean_grp(grp).mu(mu).ttt2 = All_mean_grp(grp).mu(mu).ttt1;
        
        
        
        
       
        while 1
            
            n1 = size( All_mean_grp(grp).mu(mu).ttt2( ~isnan( All_mean_grp( grp).mu( mu).ttt2( :,1) ), :), 1);
            
            sd = nanstd(All_mean_grp(grp).mu(mu).ttt2(:,1:end-1)); 
            [ysd, xsd] = max(sd);

            for suj = 1 : 12

                sujet = 1 : 12 ; 

                data = All_mean_grp(grp).mu(mu).ttt2(sujet(sujet~=suj),xsd); 

                All_mean_grp(grp).VR(mu,suj) = nanstd(data);
                All_mean_grp(grp).proba(mu,suj) = pdf('normal', All_mean_grp(grp).mu(mu).ttt2(suj,xsd), nanmean(data), nanstd(data));
            end

            seuil = 1 * 10^(-7);% mean(All_mean_grp(grp).VR(mu,:))/3; % - 3*std(Grp(grp).VR(mu,:));

            ind = All_mean_grp(grp).proba(mu,:) < seuil;

            All_mean_grp(grp).mu(mu).ttt2(sujet(ind == 1),:) = NaN;
            
            n2 = size( All_mean_grp(grp).mu(mu).ttt2( ~isnan( All_mean_grp( grp).mu( mu).ttt2( :,1) ), :), 1);
            
            if n1 == n2
                break
            end
        end
       
    end
end

%%
close all

for grp = 1 : 8
    for mu = 1 : 10
        
        figure(mu + 10*(grp-1))
        hold on
        
        data = All_mean_grp(grp).mu(mu).data(:,1:end-1);
        sd = std(data);
        
        subplot(3,1,1)
        title('BRUTE')
        
        plot(data','k')
        hold on 
        plot(mean(data), 'g', 'linewidth' , 2)
        plot(mean(data)+3*sd, 'r')
        legend(num2str(size(data,1)))
        
        %%
        
        data = All_mean_grp(grp).mu(mu).ttt1(:,1:end-1);
        data = data(~isnan(data(:,1)),:);
        sd = std(data);
        
        subplot(3,1,2)
        title('TTT1')
        
        plot(data','k')
%         plot(All_mean_grp(grp).mu(mu).data(p_suj(p_i(:)==0),1:end-1)','k')
        hold on 
        plot(mean(data), 'g', 'linewidth' , 2)
        plot(mean(data)+3*sd, 'r')
        legend(num2str(size(data,1)))
        
        %%
        data = All_mean_grp(grp).mu(mu).ttt2(:,1:end-1);
        data = data(~isnan(data(:,1)),:);
        sd = nanstd(data);
        
        subplot(3,1,3)
        title('TTT2')
        
        plot(data','k')
        hold on 
        plot(nanmean(data), 'g', 'linewidth' , 2)
        plot(nanmean(data)+3*sd, 'r')
        legend(num2str(size(data,1)))
        
        saveas(figure(mu+10*(grp-1)), ['C:\Users\p1218107\Documents\Data_Piano\graphiques\muscles\all_muscle all_grps\','grp',num2str(grp),'mu',num2str(mu), '_mean2'],'jpg')
        
    end
end

%%

keep_tot = [] ; 

for grp = 1 : 8
    for mu = 1 : 10
   
        keep_tot(mu, grp, :) = double( ~isnan( All_mean_grp(grp).mu(mu).ttt2( :,1)));
        
    end
end

