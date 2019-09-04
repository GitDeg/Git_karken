%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     Analyse statistique intragroupe                     %
%                              Protocole LA                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation

clear all;
close all;
clc;


%addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))

pathname='C:\Users\p1218107\Documents\Data_Piano\data_1.5s';

%% Ouverture des datas 

cd(pathname)

load([pathname,'\ALL_mean_1_5_manu.mat'])

%save([pathname,'\All_data_muscle.mat'],'All_data')

%% MEAN_SD cloud
close all

F=2100;


t = 1/F-525/F: 1/F : 524/(F);

for mu = 1:10
    
    if mu == 3 | mu == 7
        
        
    end
    
    if mu == 8
    end
    
    if mu == 10
    end
    
        Y1= ALL_mean(mu+1).grp1;%([3:8 , 10:12],:);
        Y2= ALL_mean(mu+1).grp2;%([3:8 , 10:12],:);
        Y3= ALL_mean(mu+1).grp3;%([3:8 , 10:12],:);
        Y4= ALL_mean(mu+1).grp4;%([3:8 , 10:12],:);
    
%     fig = figure(mu);
%     set(gcf, 'Position', get(0, 'Screensize'));
%     title(ALL_sd(mu+1).Names)
    
%     hold on
%     spm1d.plot.plot_meanSD(Y1, 'color' ,'r');
%     
%     spm1d.plot.plot_meanSD(Y2, 'color', 'g');
%     
%     spm1d.plot.plot_meanSD(Y3, 'color', 'b');
%     
%     spm1d.plot.plot_meanSD(Y4);

    
    A = [zeros(size(Y1,1),1);ones(size(Y2,1),1);zeros(size(Y3,1),1);ones(size(Y4,1),1)];
    B = [zeros(size(Y1,1),1);zeros(size(Y2,1),1);ones(size(Y3,1),1);ones(size(Y4,1),1)];
    %Subj = [[3:8 , 10:12]'; [3:8 , 10:12]'; [3:8 , 10:12]'; [3:8 , 10:12]'];
    Subj = [[1:12]' ; [1:12]' ;  [1:12]' ; [1:12]'];
    F = spm1d.stats.anova2rm([Y1; Y2; Y3; Y4],A,B,Subj);
    Fi =F.inference(0.05);
    fig = Fi.plot();
    set(gcf, 'Position', get(0, 'Screensize'));
    %Fi.plot_threshold_label();
    %Fi.plot_p_values();
    xlim([0 1049])

    
%     for i = 1 : size(Fi.clusters,2)
%         SP = round(Fi.clusters{1, i}.endpoints);
%         line([SP(1),SP(1)],ylim,'Color',[1 0 0])
%         line([SP(2),SP(2)],ylim,'Color',[1 0 0])
%     end
    
    
    %legend(ALL_sd(1).grp1,'', ALL_sd(1).grp2,'', ALL_sd(1).grp3,'', ALL_sd(1).grp4)
    
    %saveas(figure(mu),['C:\Users\p1218107\Documents\Data_Piano\graphiques\statistique12datas\mean\ANOVA2_RM\',ALL_sd(mu+1).Names,'_12anova2_RM'],'jpg')
      
end



%% Data_set

    %% Y1: Membre supérieur pressé staccato -> Grp1
    %% Y2: Membre supérieur frappé staccato -> Grp2
    %% Y3: Bassin pressé staccato -> Grp3
    %% Y4: Bassin frappé staccato -> Grp4
    %% Y5: Membre supérieur pressé tenue ->Grp5
    %% Y6: Membre supérieur frappé tenue ->Grp6
    %% Y7: Bassin pressé tenue ->Grp7
    %% Y8: Bassin frappé tenue ->Grp8

Y1= ALL_100ms(2).grp1;
Y2= ALL_100ms(2).grp2;
Y3= ALL_100ms(2).grp3;
Y4= ALL_100ms(2).grp4;
Y5= ALL_100ms(2).grp5;
Y6= ALL_100ms(2).grp6;
Y7= ALL_100ms(2).grp7;
Y8= ALL_100ms(2).grp8;


for i = 1 : 12
    sY(i,1)=size(Y1(Y1(:,size(Y1,2))==i,:),1);
    sY(i,2)=size(Y2(Y2(:,size(Y2,2))==i,:),1);
    sY(i,3)=size(Y3(Y3(:,size(Y3,2))==i,:),1);
    sY(i,4)=size(Y4(Y4(:,size(Y4,2))==i,:),1);
    sY(i,5)=size(Y5(Y5(:,size(Y5,2))==i,:),1);
    sY(i,6)=size(Y6(Y6(:,size(Y6,2))==i,:),1);
    sY(i,7)=size(Y7(Y7(:,size(Y7,2))==i,:),1);
    sY(i,8)=size(Y8(Y8(:,size(Y8,2))==i,:),1);
end

ndata=min(min(sY));
Yp1 = [];
Yp2 = [];
Yp3 = [];
Yp4 = [];
Yp5 = [];
Yp6 = [];
Yp7 = [];
Yp8 = [];

for i = 1 : 12
    Yp1 = [Yp1; datasample(Y1(Y1(:,size(Y1,2))==i,:),ndata)];
    Yp2 = [Yp2; datasample(Y2(Y2(:,size(Y2,2))==i,:),ndata)];
    Yp3 = [Yp3; datasample(Y3(Y3(:,size(Y3,2))==i,:),ndata)];
    Yp4 = [Yp4; datasample(Y4(Y4(:,size(Y4,2))==i,:),ndata)];
    Yp5 = [Yp5; datasample(Y5(Y5(:,size(Y5,2))==i,:),ndata)];
    Yp6 = [Yp6; datasample(Y6(Y6(:,size(Y6,2))==i,:),ndata)];
    Yp7 = [Yp7; datasample(Y7(Y7(:,size(Y7,2))==i,:),ndata)];
    Yp8 = [Yp8; datasample(Y8(Y8(:,size(Y8,2))==i,:),ndata)];
end

Y = [Yp1; Yp2; Yp3; Yp4; Yp5; Yp6; Yp7; Yp8];
Y_mean = [ALL_mean(2).grp1; ALL_mean(2).grp2; ALL_mean(2).grp3; ALL_mean(2).grp4];
mu = 0;

%% ONE sample t test
t = spm1d.stats.ttest(Y1, mu);

ti = t.inference(0.05,'two_tailed',false,'interp',true);    
figure(1)
subplot(2,1,1)
spm1d.plot.plot_meanSD(Y1)
subplot(2,1,2)
ti.plot()
ti.plot_threshold_label();
ti.plot_p_values();

%% Paired t test
% Incohérent car les datas ne sont pas paired
figure (2)

subplot(2,1,1)
hold on
spm1d.plot.plot_meanSD(Y1)
spm1d.plot.plot_meanSD(Y2)

subplot(2,1,2)
t = spm1d.stats.ttest_paired(Y1,Y2);
ti = t.inference(0.05,'two_tailed',false);
ti.plot()
ti.plot_threshold_label();
ti.plot_p_values();

%% Two-sample t test
close all
iterations = 500; 

for mu = 1 : 10
    
    Y1= ALL_mean(mu+1).grp1([3:8 , 10:12],:);
    Y2= ALL_mean(mu+1).grp2([3:8 , 10:12],:);
    Y3= ALL_mean(mu+1).grp3([3:8 , 10:12],:);
    Y4= ALL_mean(mu+1).grp4([3:8 , 10:12],:);
    Subj = [3:8 , 10:12]';%;[3:8 , 10:12]';[3:8 , 10:12]';[3:8 , 10:12]'];
    
    t12 = spm1d.stats.nonparam.ttest_paired(Y1,Y2); %,'equal_var',false
    t12i = t12.inference(0.05,'two_tailed',true, 'iterations', iterations);
    t13 = spm1d.stats.nonparam.ttest_paired(Y1,Y3);
    t13i = t13.inference(0.05,'two_tailed',true, 'iterations', iterations);

    t24 = spm1d.stats.nonparam.ttest_paired(Y2,Y4);
    t24i = t24.inference(0.05,'two_tailed',true, 'iterations', iterations);
    
    t34 = spm1d.stats.nonparam.ttest_paired(Y3,Y4);
    t34i = t34.inference(0.05,'two_tailed',true, 'iterations', iterations);
    
    t43 = spm1d.stats.nonparam.ttest_paired(Y4,Y3);
    t43i = t43.inference(0.05,'two_tailed',true, 'iterations', iterations);
    
    %% Plot 1
    fig = figure (mu);
    set(gcf, 'Position', get(0, 'Screensize'));
    subplot(3,1,1)
    title('MS : Frappé vs pressé')
    hold on
    spm1d.plot.plot_meanSD(Y1, 'color' ,'r');
    spm1d.plot.plot_meanSD(Y2, 'color', 'g');
    t12i.plot()
    for i = 1 : size(t12i.clusters,2)
        SP = round(t12i.clusters{1, i}.endpoints);
        line([SP(1),SP(1)],ylim,'Color',[1 0 0])
        line([SP(2),SP(2)],ylim,'Color',[1 0 0])
    end
    legend('Pressé','', 'Frappé','', 'p value')

    %%
    subplot(3,2,5)
    title('Frappé : Bassin vs MS')
    hold on
    spm1d.plot.plot_meanSD(Y1, 'color' ,'r');
    spm1d.plot.plot_meanSD(Y3, 'color', 'b');
    t13i.plot()
    for i = 1 : size(t13i.clusters,2)
        SP = round(t13i.clusters{1, i}.endpoints);
        line([SP(1),SP(1)],ylim,'Color',[1 0 0])
        line([SP(2),SP(2)],ylim,'Color',[1 0 0])
    end
   legend('MS','', 'Bassin','', 'p value')
   
  %%
    subplot(3,2,6)
    title('Pressé : Bassin vs MS')
    hold on
    spm1d.plot.plot_meanSD(Y2, 'color' ,'g');
    spm1d.plot.plot_meanSD(Y4, 'color', 'k');
    t24i.plot()
    for i = 1 : size(t24i.clusters,2)
        SP = round(t24i.clusters{1, i}.endpoints);
        line([SP(1),SP(1)],ylim,'Color',[1 0 0])
        line([SP(2),SP(2)],ylim,'Color',[1 0 0])
    end
    legend('MS','', 'Bassin','', 'p value')

    
    %%     
    subplot(3,1,2)
    title('Bassin : Frappé vs pressé')
    hold on
    spm1d.plot.plot_meanSD(Y3, 'color' ,'b');
    spm1d.plot.plot_meanSD(Y4, 'color', 'k');
    t34i.plot()
    for i = 1 : size(t34i.clusters,2)
        SP = round(t34i.clusters{1, i}.endpoints);
        line([SP(1),SP(1)],ylim,'Color',[1 0 0])
        line([SP(2),SP(2)],ylim,'Color',[1 0 0])
    end
    legend('Pressé','', 'Frappé','', 'p value')

    saveas(fig,['C:\Users\p1218107\Documents\Data_Piano\graphiques\statistique12datas\mean\ttest_paired\',ALL_mean(mu+1).Names,'_ttest_paired'],'jpg')
    
end

%% One-way ANOVA 
figure(4)

Y1= ALL_(2).grp1;
Y2= ALL_mean(2).grp2;
Y3= ALL_mean(2).grp3;
Y4= ALL_mean(2).grp4;

A = [zeros(size(Y1,1),1);ones(size(Y2,1),1);2*ones(size(Y3,1),1);3*ones(size(Y4,1),1)];

F = spm1d.stats.nonparam.anova1rm([Y1; Y2; Y3; Y4],A);
Fi =F.inference(0.05,'iterations', 200);
Fi.plot()
Fi.plot_threshold_label();
Fi.plot_p_values();

%% One-way ANOVA repeated-measure 

figure(4)

A = [zeros(size(Y1,1),1); ones(size(Y2,1),1); 2*ones(size(Y3,1),1); 3*ones(size(Y4,1),1)];%...
     %4*ones(size(Yp2,1),1);5*ones(size(Yp3,1),1);6*ones(size(Yp4,1),1); 7*ones(size(Yp4,1),1)];

iterations = 100;

F  = spm1d.stats.nonparam.anova1rm(Y(:,1:size(Y,2)-1), A, Subj); %Y(:,size(Y,2))
Fi = F.inference(0.05,'iterations', iterations);
Fi.plot()
Fi.plot_threshold_label();
Fi.plot_p_values();

%% Two-way ANOVA


for mu = 1 : 10

    Y1 = ALL_mean_ali_1_5(mu + 1).grp1(:,1:end-1); 
    Y2 = ALL_mean_ali_1_5(mu + 1).grp2(:,1:end-1);
    Y3 = ALL_mean_ali_1_5(mu + 1).grp3(:,1:end-1);
    Y4 = ALL_mean_ali_1_5(mu + 1).grp4(:,1:end-1);

    A=[zeros(12,1);ones(12,1);zeros(12,1);ones(12,1)];
    Subj=[[1:12]'; [1:12]'; [1:12]'; [1:12]'];

    
    alpha      = 0.05;
    iterations = 500;
    F = spm1d.stats.nonparam.anova2([Y1;Y2;Y3;Y4],A,Subj);
    Fi = F.inference(alpha,'iterations',iterations);
    
    fig = Fi.plot();
    title(subplot(2,2,1),'Main Pressé / Frappé')
    title(subplot(2,2,2),'Main Bassin / Membre Supérieur')
    title(subplot(2,2,3),'Interaction')
    subplot(2,2,4)
    hold on
    spm1d.plot.plot_meanSD(Y1, 'color' ,'r');
    spm1d.plot.plot_meanSD(Y2, 'color', 'g');
    spm1d.plot.plot_meanSD(Y3, 'color', 'b');
    spm1d.plot.plot_meanSD(Y4);
    title('Courbes Mean SD')
    
    set(gcf, 'Position', get(0, 'Screensize'))
    suptitle(ALL_mean_1_5_manu(mu+1).Names) 

    legend('MS pressé','', 'MS frappé','', 'Bassin pressé','', 'Bassin frappé')
    saveas(figure(mu),['C:\Users\p1218107\Documents\Data_Piano\graphiques\statistique\EMGmoyen\ANOVA2_non_para_data_manu\',ALL_mean_1_5_manu(mu+1).Names,'_anova2_data_manu'],'jpg')

end

%% Hotelling Data_set

ALL_mean(9).grp2(12,:) = mean(ALL_mean(9).grp2);
ALL_mean(11).grp4(12,:) = mean(ALL_mean(11).grp4);
ALL_mean(8).grp4(12,:) = mean(ALL_mean(8).grp4);
ALL_mean(4).grp4(12,:) = mean(ALL_mean(4).grp4);

GRP1 = [];
GRP2 = [];
GRP3 = [];
GRP4 = [];


for i = 1:10
    GRP1(:,:,i) = ALL_mean(i+1).grp1;
    GRP2(:,:,i) = ALL_mean(i+1).grp2;
    GRP3(:,:,i) = ALL_mean(i+1).grp3;
    GRP4(:,:,i) = ALL_mean(i+1).grp4;
end

%% One-sample Hotelling's T² test
figure(7)

T2 = spm1d.stats.hotellings(GRP1(:,:,[2,3,5,6]), mu);
T2i = T2.inference(0.05);
T2i.plot()
T2i.plot_threshold_label();
T2i.plot_p_values();

%% Paired Hotelling's T² test

figure(8)
T2 = spm1d.stats.hotellings_paired(GRP1(:,:,[2,3,5,6]),GRP3(:,:,[2,3,5,6]));
T2i = T2.inference(0.05);
T2i.plot()
T2i.plot_threshold_label();
T2i.plot_p_values();

%% Two-sample Hotelling's T² test

figure(8)
T2 = spm1d.stats.hotellings2(cat(3,GRP1(:,:,[2,3,5,6]),GRP2(:,:,[2,3,5,6])),cat(3,GRP3(:,:,[2,3,5,6]),GRP4(:,:,[2,3,5,6])));
T2i = T2.inference(0.05);
T2i.plot()
T2i.plot_threshold_label();
T2i.plot_p_values();

%%

