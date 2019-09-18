%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                          ANOVA 2 RM                               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc 

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
%%
% pathname = '\\10.89.24.15\j\Valentin\SAUT';
pathname = 'C:\Users\p1218107\Documents\Data_Piano\SAUT';

sujet = {'\001','\002','\003','\004','\005','\006','\007','\008','\009','\010','\011','\012'};

load([pathname '\All_data.mat'])

%%
close all

% F=2100;
% t = 1/F -1000/(2*F) : 1/F : 1000/(2*F);

ordre_mu = [9 8 2 1 3 4 7 5 6];

for mu = ordre_mu
       
    suj_set = [1:12]';
    
    Y1 = All_data(1).EMG(1).mean_EMG(mu).data_MVC;
    Y2 = All_data(2).EMG(1).mean_EMG(mu).data_MVC;
    Y3 = All_data(3).EMG(1).mean_EMG(mu).data_MVC;
    Y4 = All_data(4).EMG(1).mean_EMG(mu).data_MVC;
    
    suj_set = suj_set(~isnan(Y1(suj_set,1)));
    suj_set = suj_set(~isnan(Y2(suj_set,1)));
    suj_set = suj_set(~isnan(Y3(suj_set,1)));
    suj_set = suj_set(~isnan(Y4(suj_set,1)));
    

    Y = [Y1(suj_set,:); Y2(suj_set,:); Y3(suj_set,:); Y4(suj_set,:)];
    
    iterations = 1000;
    
    A = [zeros(size(Y1(suj_set,:),1),1);ones(size(Y2(suj_set,:),1),1);zeros(size(Y3(suj_set,:),1),1);ones(size(Y4(suj_set,:),1),1)];  % 0 = sans bassin   / 1 = avec bassin
    B = [zeros(size(Y1(suj_set,:),1),1);zeros(size(Y2(suj_set,:),1),1);ones(size(Y3(suj_set,:),1),1);ones(size(Y4(suj_set,:),1),1)];  % 0 = anti-horaire  / 1 = Demi-cercle
    
    Subj = repmat(suj_set,4,1);
    
    F = spm1d.stats.nonparam.anova2rm(Y,A,B,Subj);
    Fi =F.inference(0.05,'iteration',iterations);
    
    fig = Fi.plot();

    title(subplot(2,2,1),'Bassin implication')
    title(subplot(2,2,2),'anti-horaire / demi-cercle')
    title(subplot(2,2,3),'Interaction')
    subplot(2,2,4)
    hold on
    spm1d.plot.plot_meanSD(1:length(Y1),Y1(suj_set,:), 'color', 'r');
    spm1d.plot.plot_meanSD(1:length(Y2),Y2(suj_set,:), 'color', 'g');
    spm1d.plot.plot_meanSD(1:length(Y3),Y3(suj_set,:), 'color', 'b');
    spm1d.plot.plot_meanSD(1:length(Y4),Y4(suj_set,:));
    title('Courbes Mean SD')
    
%     set(gcf, 'Position', get(0, 'Screensize'))
    suptitle(All_data(1).EMG(1).mu(mu).name) 

    legend('Anti','', 'Anti-bassin','', 'Demi','', 'Demi-bassin')   
end

%%