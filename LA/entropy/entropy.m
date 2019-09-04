%% entropy - synergy analysis

clear all;
close all;
clc;

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))

%% entropy of each muscle


for suj = 1 : 12
    for mui = 1 : 10
       for muj = 1 : 10

            H.grp1(suj,mui,muj) = mutualinfo(ALL_mean_1_5_manu(mui+1).grp1(suj,1:end-1), ALL_mean_1_5_manu(muj+1).grp1(suj,1:end-1))/...
                jointentropy(ALL_mean_1_5_manu(mui+1).grp1(suj,1:end-1), ALL_mean_1_5_manu(muj+1).grp1(suj,1:end-1));
            
            H.grp2(suj,mui,muj) = mutualinfo(ALL_mean_1_5_manu(mui+1).grp2(suj,1:end-1), ALL_mean_1_5_manu(muj+1).grp2(suj,1:end-1))/...
                jointentropy(ALL_mean_1_5_manu(mui+1).grp2(suj,1:end-1), ALL_mean_1_5_manu(muj+1).grp2(suj,1:end-1));
            
            H.grp3(suj,mui,muj) = mutualinfo(ALL_mean_1_5_manu(mui+1).grp3(suj,1:end-1), ALL_mean_1_5_manu(muj+1).grp3(suj,1:end-1))/...
                jointentropy(ALL_mean_1_5_manu(mui+1).grp3(suj,1:end-1), ALL_mean_1_5_manu(muj+1).grp3(suj,1:end-1));
            
            H.grp4(suj,mui,muj) = mutualinfo(ALL_mean_1_5_manu(mui+1).grp4(suj,1:end-1), ALL_mean_1_5_manu(muj+1).grp4(suj,1:end-1))/...
                jointentropy(ALL_mean_1_5_manu(mui+1).grp4(suj,1:end-1), ALL_mean_1_5_manu(muj+1).grp4(suj,1:end-1));
       
       end
    end
end



%%

for mu = 1 : 10
    subplot(5,2,mu)
    pos = [(1:9)-0.15, (1:9)-0.10, (1:9)+0.10, (1:9)+0.15];
    col = [repmat({'r'},1,9) repmat({'g'},1,9) repmat({'b'},1,9) repmat({'k'},1,9)];
    
    entro_mu = 1:10;
    entro_name = {'Flé_1';'Flé_2';'Ext';'Bcps';'Tcps';'DltAnt';'DltMed';'GrdPec';'TrpSup';'GrdDent'} ; 
    
    boxplot([H.grp1(:,entro_mu(entro_mu~=mu),mu) H.grp2(:,entro_mu(entro_mu~=mu),mu) H.grp3(:,entro_mu(entro_mu~=mu),mu) H.grp4(:,entro_mu(entro_mu~=mu),mu)],...
        'position', pos ,'PlotStyle','compact', 'ColorGroup', col, 'LabelOrientation', 'inline')
    ylim([0 0.6])
    ylabel(ALL_mean_1_5_manu(mu+1).Names)
    set(gca, 'xtick', [1 2 3 4 5 6 7 8 9])
    set(gca, 'xticklabel', entro_name(entro_mu(entro_mu~=mu)))

    %groups={[3,6],[1,2],[5,4]};
    %sigstar(groups,[0.001,0.05,0.04]);
    
end
