clear all;
close all;
clc;


%addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB\spm1dmatlab-master'))

pathname='C:\Users\p1218107\Documents\Data_Piano';
%% Ouverture des datas 

cd(pathname)

load([pathname,'\All_data_muscle_30ms.mat'])

%save([pathname,'\All_data_muscle_30ms.mat'],'ALL_30ms')

%% ANOVA 


for mu = [1 2 3 5 9] 
    
    load(['C:\Users\p1218107\Documents\Data_Piano\muscle_interet\',ALL_30ms(mu+1).Names,'_F.mat'])
    load(['C:\Users\p1218107\Documents\Data_Piano\muscle_interet\',ALL_30ms(mu+1).Names,'_Fi.mat'])
    fig = figure(mu);
    set(gcf, 'Position', get(0, 'Screensize'));
    title(ALL_30ms(mu+1).Names)
    subplot(2,1,1)
    hold on
    spm1d.plot.plot_meanSD(ALL_30ms(mu+1).grp1(:,1:end-1), 'color' ,'r');
    
    spm1d.plot.plot_meanSD(ALL_30ms(mu+1).grp2(:,1:end-1), 'color', 'g');
    
    spm1d.plot.plot_meanSD(ALL_30ms(mu+1).grp3(:,1:end-1), 'color', 'b');
    
    spm1d.plot.plot_meanSD(ALL_30ms(mu+1).grp4(:,1:end-1));
    
    for i = 1 : size(Fi.clusters,2)
        bar(round(Fi.clusters{1, i}.endpoints),[8,8],0.01)
    end
    
    
    legend(ALL_30ms(1).grp1,'', ALL_30ms(1).grp2,'', ALL_30ms(1).grp3,'', ALL_30ms(1).grp4)
    xlim([0 1049])

    subplot(2,1,2)
    Fi.plot()
    xlim([0 1049])
    Fi.plot_threshold_label();
    Fi.plot_p_values();
    
     

    %saveas(fig,['C:\Users\p1218107\Documents\Data_Piano\graphiques\ANOVA RM 30ms\interet\',ALL_30ms(mu+1).Names,'clust'],'jpeg')
end
%%

