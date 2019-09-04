%% Initialisation

clear all;
close all;
clc;


%addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB\spm1dmatlab-master'))

pathname='C:\Users\p1218107\Documents\Data_Piano';

%% Ouverture des datas 


cd(pathname)

load([pathname,'\All_data_muscle_30ms.mat'])

%%
close all

F=2100;

sujets = 1 : 12;
ndata = [];
t = 1/F-525/F: 1/F : 524/(F);

for mu = [1, 2, 4:6, 9]
    for i = sujets
        
        sY(i,1)=size(ALL_30ms(mu+1).grp1(ALL_30ms(mu+1).grp1(:,size(Y1,2))==i,:),1);
        sY(i,2)=size(ALL_30ms(mu+1).grp2(ALL_30ms(mu+1).grp2(:,size(Y2,2))==i,:),1);
        sY(i,3)=size(ALL_30ms(mu+1).grp3(ALL_30ms(mu+1).grp3(:,size(Y3,2))==i,:),1);
        sY(i,4)=size(ALL_30ms(mu+1).grp4(ALL_30ms(mu+1).grp4(:,size(Y4,2))==i,:),1);

    end
    ndata=[ndata , nanmin(nanmin(sY))];    
end

%%

Yp1 = [];
Yp2 = [];
Yp3 = [];
Yp4 = [];

for mu = [1, 2, 4:6, 9]
    
    Y1 = [];
    Y2 = [];
    Y3 = [];
    Y4 = [];

    
    Y1= ALL_30ms(mu+1).grp1;
    Y2= ALL_30ms(mu+1).grp2;
    Y3= ALL_30ms(mu+1).grp3;
    Y4= ALL_30ms(mu+1).grp4;
   
    
    for i = sujets
    Yp1 = [Yp1; datasample(Y1(Y1(:,size(Y1,2))==i,:),min(ndata))];
    Yp2 = [Yp2; datasample(Y2(Y2(:,size(Y2,2))==i,:),min(ndata))];
    Yp3 = [Yp3; datasample(Y3(Y3(:,size(Y3,2))==i,:),min(ndata))];
    Yp4 = [Yp4; datasample(Y4(Y4(:,size(Y4,2))==i,:),min(ndata))];
    end
end

%%
    Y = [Yp1; Yp2; Yp3; Yp4];

    iterations = 20 ; %1000;
    A = [zeros(size(Yp1,1),1);ones(size(Yp2,1),1);zeros(size(Yp3,1),1);ones(size(Yp4,1),1)];
    B = [zeros(size(Yp1,1),1);zeros(size(Yp2,1),1);ones(size(Yp3,1),1);ones(size(Yp4,1),1)];
    C = [0*ones(min(ndata),1); 1*ones(min(ndata),1); 2*ones(min(ndata),1);...
         3*ones(min(ndata),1); 4*ones(min(ndata),1); 5*ones(min(ndata),1);];
    C = [C; C; C; C; C; C; C; C; C; C; C; C];
    C = [C; C; C; C];
    
    Subj = Y(:,size(Y,2));
    
    F = spm1d.stats.anova3rm(Y(:,1:1049),A,B,C,Subj);
    
    Fi =F.inference(0.05,'iteration',iterations);
    
    fig = Fi.plot();

    title(subplot(2,2,1),'Main Pressé / Frappé')
    title(subplot(2,2,2),'Main Bassin / Membre Supérieur')
    title(subplot(2,2,3),'Interaction')
    subplot(2,2,4)
    hold on
    spm1d.plot.plot_meanSD(Y1(:,1:1049), 'color' ,'r');
    spm1d.plot.plot_meanSD(Y2(:,1:1049), 'color', 'g');
    spm1d.plot.plot_meanSD(Y3(:,1:1049), 'color', 'b');
    spm1d.plot.plot_meanSD(Y4(:,1:1049));
    title('Courbes Mean SD')
    
    set(gcf, 'Position', get(0, 'Screensize'))
    suptitle(ALL_30ms(mu+1).Names) 

    legend('MS pressé','', 'MS frappé','', 'Bassin pressé','', 'Bassin frappé')

    xlim([0 1049])

%     save(['C:\Users\p1218107\Documents\Data_Piano\graphiques\ANOVA2 RM all_data\muscles\',ALL_30ms(mu+1).Names,'_anova2_RM_F'],'F')
%     save(['C:\Users\p1218107\Documents\Data_Piano\graphiques\ANOVA2 RM all_data\muscles\',ALL_30ms(mu+1).Names,'_anova2_RM_Fi'],'Fi')
%     saveas(figure(mu),['C:\Users\p1218107\Documents\Data_Piano\graphiques\ANOVA2 RM all_data\',ALL_30ms(mu+1).Names,'_anova2_RM_alldata'],'jpg')
%       


