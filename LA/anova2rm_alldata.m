%% Initialisation

clear all;
close all;
clc;


addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))

pathname='C:\Users\p1218107\Documents\Data_Piano\data_1.5s';

%% Ouverture des datas 

cd(pathname)

load([pathname,'\ALL_data_1_5_manu.mat'])

%%
close all

F=2100;


t = 1/F-525/F: 1/F : 524/(F);

for mu = 2%:10
    
    sujets = [];
    clear('sY','Y','Y1','Y2','Y3','Y4','Yp1','Yp2','Yp3','Yp4','A','B')
    
    if mu == 3 | mu == 7
        Y1= ALL_data_1_5_manu(mu+1).grp1(ALL_data_1_5_manu(mu+1).grp1(:,size(ALL_data_1_5_manu(mu+1).grp1,2))~=9,:);
        Y2= ALL_data_1_5_manu(mu+1).grp2(ALL_data_1_5_manu(mu+1).grp2(:,size(ALL_data_1_5_manu(mu+1).grp1,2))~=9,:);
        Y3= ALL_data_1_5_manu(mu+1).grp3(ALL_data_1_5_manu(mu+1).grp3(:,size(ALL_data_1_5_manu(mu+1).grp1,2))~=9,:);
        Y4= ALL_data_1_5_manu(mu+1).grp4(ALL_data_1_5_manu(mu+1).grp4(:,size(ALL_data_1_5_manu(mu+1).grp1,2))~=9,:);
        sujets = [1:8 , 10:12];
        
    elseif mu == 8
        Y1= ALL_data_1_5_manu(mu+1).grp1(ALL_data_1_5_manu(mu+1).grp1(:,size(ALL_data_1_5_manu(mu+1).grp1,2))~=2,:);
        Y2= ALL_data_1_5_manu(mu+1).grp2(ALL_data_1_5_manu(mu+1).grp2(:,size(ALL_data_1_5_manu(mu+1).grp1,2))~=2,:);
        Y3= ALL_data_1_5_manu(mu+1).grp3(ALL_data_1_5_manu(mu+1).grp3(:,size(ALL_data_1_5_manu(mu+1).grp1,2))~=2,:);
        Y4= ALL_data_1_5_manu(mu+1).grp4(ALL_data_1_5_manu(mu+1).grp4(:,size(ALL_data_1_5_manu(mu+1).grp1,2))~=2,:);
        sujets = [1, 3:12];
        
    elseif mu == 10
        Y1= ALL_data_1_5_manu(mu+1).grp1(ALL_data_1_5_manu(mu+1).grp1(:,size(ALL_data_1_5_manu(mu+1).grp1,2))~=1,:);
        Y2= ALL_data_1_5_manu(mu+1).grp2(ALL_data_1_5_manu(mu+1).grp2(:,size(ALL_data_1_5_manu(mu+1).grp1,2))~=1,:);
        Y3= ALL_data_1_5_manu(mu+1).grp3(ALL_data_1_5_manu(mu+1).grp3(:,size(ALL_data_1_5_manu(mu+1).grp1,2))~=1,:);
        Y4= ALL_data_1_5_manu(mu+1).grp4(ALL_data_1_5_manu(mu+1).grp4(:,size(ALL_data_1_5_manu(mu+1).grp1,2))~=1,:);
        sujets = 2:12;
    else
    
        Y1= ALL_data_1_5_manu(mu+1).grp1;
        Y2= ALL_data_1_5_manu(mu+1).grp2;
        Y3= ALL_data_1_5_manu(mu+1).grp3;
        Y4= ALL_data_1_5_manu(mu+1).grp4;
        sujets = 1 : 12;
    end
    
    sY= nan(12,4);
    for i = sujets
        sY(i,1)=size(Y1(Y1(:,size(Y1,2))==i,:),1);
        sY(i,2)=size(Y2(Y2(:,size(Y2,2))==i,:),1);
        sY(i,3)=size(Y3(Y3(:,size(Y3,2))==i,:),1);
        sY(i,4)=size(Y4(Y4(:,size(Y4,2))==i,:),1);
    end

    ndata=nanmin(nanmin(sY));
    Yp1 = [];
    Yp2 = [];
    Yp3 = [];
    Yp4 = [];

    for i = sujets
        Yp1 = [Yp1; datasample(Y1(Y1(:,size(Y1,2))==i,:),ndata)];
        Yp2 = [Yp2; datasample(Y2(Y2(:,size(Y2,2))==i,:),ndata)];
        Yp3 = [Yp3; datasample(Y3(Y3(:,size(Y3,2))==i,:),ndata)];
        Yp4 = [Yp4; datasample(Y4(Y4(:,size(Y4,2))==i,:),ndata)];
    end

    Y = [Yp1; Yp2; Yp3; Yp4];

    iterations = 1000;
    A = [zeros(size(Yp1,1),1);ones(size(Yp2,1),1);zeros(size(Yp3,1),1);ones(size(Yp4,1),1)];
    B = [zeros(size(Yp1,1),1);zeros(size(Yp2,1),1);ones(size(Yp3,1),1);ones(size(Yp4,1),1)];
    Subj = Y(:,size(Y,2));
    
    F = spm1d.stats.nonparam.anova2rm(Y(:,1:size(Y,2)-1),A,B,Subj);
    Fi =F.inference(0.05,'iteration',iterations);
    
    fig = Fi.plot();

    title(subplot(2,2,1),'Main Pressé / Frappé')
    title(subplot(2,2,2),'Main Bassin / Membre Supérieur')
    title(subplot(2,2,3),'Interaction')
    subplot(2,2,4)
    hold on
    spm1d.plot.plot_meanSD(Y1(:,1:size(Y,2)-1), 'color' ,'r');
    spm1d.plot.plot_meanSD(Y2(:,1:size(Y,2)-1), 'color', 'g');
    spm1d.plot.plot_meanSD(Y3(:,1:size(Y,2)-1), 'color', 'b');
    spm1d.plot.plot_meanSD(Y4(:,1:size(Y,2)-1));
    title('Courbes Mean SD')
    
    set(gcf, 'Position', get(0, 'Screensize'))
    suptitle(ALL_data_1_5_manu(mu+1).Names) 

    legend('MS pressé','', 'MS frappé','', 'Bassin pressé','', 'Bassin frappé')

    xlim([0 size(Y,2)-1])

    save(['C:\Users\p1218107\Documents\Data_Piano\graphiques\ANOVA2 RM all_data\1.5s_manu\muscles\',ALL_data_1_5_manu(mu+1).Names,'_anova2_RM_F'],'F')
    save(['C:\Users\p1218107\Documents\Data_Piano\graphiques\ANOVA2 RM all_data\1.5s_manu\muscles\',ALL_data_1_5_manu(mu+1).Names,'_anova2_RM_Fi'],'Fi')
    saveas(figure(1),['C:\Users\p1218107\Documents\Data_Piano\graphiques\ANOVA2 RM all_data\1.5s_manu\',ALL_data_1_5_manu(mu+1).Names,'_anova2_RM_alldata'],'jpg')
%       
end

