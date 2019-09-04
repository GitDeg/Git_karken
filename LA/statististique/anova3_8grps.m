%% Initialisation

clear all;
close all;
clc;


%addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))

pathname='C:\Users\p1218107\Documents\Data_Piano\data_newREPS';

%% Ouverture des datas 


cd(pathname)

%load([pathname,'\All_grp.mat'])
load([pathname,'\All_data2.mat'])


ordre_mu = [9 8 2 1 3 4 7 5 6];

ordre_mu_name = {'Triceps'; 'Biceps'; 'Anterior Deltoid'; 'Middle Deltoid'; 'Upper Trapezius'; 'Serratus Anterior'; 'Great Pectoral'; 'Extensor Digitorum'; 'Flexor Digitorum' };

condition_name = {'All-body Struck Staccato'; 'All-body Struck Tenuto'; 'All-body Pressed Staccato'; 'All-body Pressed Tenuto';...
                  'Shoulder Struck Staccato'; 'Shoulder Struck Tenuto'; 'Shoulder Pressed Staccato';'Shoulder Pressed Tenuto'};


%%
close all

F=2100/5;
t = 1/F -630/(2*F) : 1/F : 630/(2*F);
%t = 1 : 630 ;
% F=150;
% 
% t = 1/F -227/(2*F) : 1/F : 227/(2*F);

for mu = ordre_mu(3) %1% : 9
    
    sujets = [];
    clear('sY','Y','Y1','Y2','Y3','Y4','Yp1','Yp2','Yp3','Yp4','A','B')
    for grp = 1 : 8
        suff = num2str(grp);
        
%         assignin('base', sprintf('Y%s',suff), All_data(grp).suj_mean(mu).mu_ssAberrant)
%         assignin('base', sprintf('Y%s',suff), All_data(grp).suj_mean(mu).mu_maj)
         assignin('base', sprintf('Y%s',suff), [All_data(grp).suj_mean(mu).mu, [1:12]'])
%            assignin('base', sprintf('Y%s',suff), Cine_touche(grp).grp)
    end
    
    Y1(isnan(Y1(:,1)),:) = [];
    Y2(isnan(Y2(:,1)),:) = [];
    Y3(isnan(Y3(:,1)),:) = [];
    Y4(isnan(Y4(:,1)),:) = [];
    Y5(isnan(Y5(:,1)),:) = [];
    Y6(isnan(Y6(:,1)),:) = [];
    Y7(isnan(Y7(:,1)),:) = [];
    Y8(isnan(Y8(:,1)),:) = [];
    
    suj_set = mintersect(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8);

    Y1 = Y1(ismember(Y1(:,end),suj_set),1:end-1);
    Y2 = Y2(ismember(Y2(:,end),suj_set),1:end-1);
    Y3 = Y3(ismember(Y3(:,end),suj_set),1:end-1);
    Y4 = Y4(ismember(Y4(:,end),suj_set),1:end-1);
    Y5 = Y5(ismember(Y5(:,end),suj_set),1:end-1);
    Y6 = Y6(ismember(Y6(:,end),suj_set),1:end-1);
    Y7 = Y7(ismember(Y7(:,end),suj_set),1:end-1);
    Y8 = Y8(ismember(Y8(:,end),suj_set),1:end-1);
   
    
    Y = [Y1; Y2; Y3; Y4; Y5; Y6; Y7; Y8];

    iterations = 1000 ; %1000;
    A = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); ones(size(Y3,1),1); ones(size(Y4,1),1); ...
        zeros(size(Y5,1),1); zeros(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % frappé (0) / pressé (1)
    
    B = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); zeros(size(Y3,1),1); zeros(size(Y4,1),1); ...
        ones(size(Y5,1),1); ones(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % Bassin (0) / MS (1)
    
    C = [zeros(size(Y1,1),1); ones(size(Y2,1),1); zeros(size(Y3,1),1); ones(size(Y4,1),1); ...
        zeros(size(Y5,1),1); ones(size(Y6,1),1); zeros(size(Y7,1),1); ones(size(Y8,1),1)]; % Sta (0) / Tenue (1)
    
    Subj = [suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set]; %Y(:,size(Y,2));

    F = spm1d.stats.nonparam.anova3rm(Y, A, B, C,Subj);
    Fi =F.inference(0.05,'iteration',iterations);

    
%% preparation data ANOVA3 non répété
%     
%     
%     
%     Y1 = Y1(randperm(size(Y1,1),7)', 1:end-1);
%     Y2 = Y2(randperm(size(Y2,1),7)', 1:end-1);
%     Y3 = Y3(randperm(size(Y3,1),7)', 1:end-1);
%     Y4 = Y4(randperm(size(Y4,1),7)', 1:end-1);
%     Y5 = Y5(randperm(size(Y5,1),7)', 1:end-1);
%     Y6 = Y6(randperm(size(Y6,1),7)', 1:end-1);
%     Y7 = Y7(randperm(size(Y7,1),7)', 1:end-1);
%     Y8 = Y8(randperm(size(Y8,1),7)', 1:end-1);
%    
%     
%     Y = [Y1; Y2; Y3; Y4; Y5; Y6; Y7; Y8];
% 
%     iterations = 1000 ; %1000;
%     A = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); ones(size(Y3,1),1); ones(size(Y4,1),1); ...
%         zeros(size(Y5,1),1); zeros(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % frappé (0) / pressé (1)
%     
%     B = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); zeros(size(Y3,1),1); zeros(size(Y4,1),1); ...
%         ones(size(Y5,1),1); ones(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % Bassin (0) / MS (1)
%     
%     C = [zeros(size(Y1,1),1); ones(size(Y2,1),1); zeros(size(Y3,1),1); ones(size(Y4,1),1); ...
%         zeros(size(Y5,1),1); ones(size(Y6,1),1); zeros(size(Y7,1),1); ones(size(Y8,1),1)]; % Sta (0) / Tenue (1)
%     
%     F = spm1d.stats.nonparam.anova3(Y, A, B, C);
%     Fi =F.inference(0.05,'iteration',iterations);

    
    %%
    figure(mu)
    %set(gcf, 'Position', get(0, 'Screensize'))
    subplot(3,3,8:9)
    title('Means All grps')
    hold on
    plot(t, mean(Y1), 'r', 'LineWidth', 1.5)
    plot(t, mean(Y2), 'g', 'LineWidth', 1.5)
    plot(t, mean(Y3), 'b', 'LineWidth', 1.5)
    plot(t, mean(Y4), 'k', 'LineWidth', 1.5)
    plot(t, mean(Y5), ':r', 'LineWidth', 1.5)
    plot(t, mean(Y6), ':g', 'LineWidth', 1.5)
    plot(t, mean(Y7), ':b', 'LineWidth', 1.5)
    plot(t, mean(Y8), ':k', 'LineWidth', 1.5)
    
    
    Titles = {'A. Struck/Pressed' 'B. All-body/Shoulder' 'C. Staccato/Tenuto' 'Interaction A/B' 'Interaction A/C' 'Interaction B/C' 'Interaction A/B/C'};
    
    for clu = 1 : 7
        subplot(3,3,clu)
        plot(t, Fi.SPMs{1, clu}.z,'k')
        hold on 
        line(xlim , [Fi.SPMs{1, clu}.zstar Fi.SPMs{1, clu}.zstar] ,'Color',[1 0 0], 'LineStyle','--')
        title(subplot(3,3,clu), Titles(clu))
        
        for i = 1 : size(Fi.SPMs{1, clu}.clusters,2)
            SP = round(Fi.SPMs{1, clu}.clusters{1, i}.endpoints);
            if SP(1) == 0
                SP(1) = 1;
            end
            subplot(3,3,8:9)
            hold on
        	line(t([SP(1),SP(1)]),ylim,'Color',[1 0 0])
            line(t([SP(2),SP(2)]),ylim,'Color',[1 0 0])
        end
    end
    
 %  suptitle(ordre_mu_name(mu))
   
 %  saveas(figure(mu), ['C:\Users\p1218107\Documents\Data_Piano\graphiques\New_reps\Anova3\10 cycles\' All_data(grp).suj_mean(mu).name '_all_interactions'], 'jpeg')
    
end

