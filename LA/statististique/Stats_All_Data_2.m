%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     Analyse statistique intragroupe                     %
%                              Protocole LA                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation

clear all;
close all;
clc;


%addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB\spm1dmatlab-master'))

pathname='C:\Users\p1218107\Documents\Data_Piano\data_1.5s';

%% Ouverture des datas 

cd(pathname)

load([pathname,'\ALL_mean_1_5_manu.mat'])

%save([pathname,'\All_data_muscle.mat'],'All_data')

%% Dataset

mu=11;
Y1= ALL_30ms(mu).grp1(ALL_30ms(mu).grp1(:,size(ALL_30ms(mu).grp1,2))==1, 1: size(ALL_30ms(mu).grp1,2)-1 ) ;
Y2= ALL_30ms(mu).grp1(ALL_30ms(mu).grp1(:,size(ALL_30ms(mu).grp1,2))==3, 1: size(ALL_30ms(mu).grp1,2)-1 ) ;


%% TTEST
figure

subplot(2,1,1)
hold on
spm1d.plot.plot_meanSD(Y1)
spm1d.plot.plot_meanSD(Y2)

subplot(2,1,2)
%t = spm1d.stats.ttest2([Y1;Y3],[Y2;Y4]);
t = spm1d.stats.ttest2(Y1,Y2);
ti = t.inference(0.05,'two_tailed',true, 'interp',true);
ti.plot()
ti.plot_threshold_label();
ti.plot_p_values();

%% ANOVA

mu= 5;
figure
for suj = 1 : 12

    Y1= ALL_30ms(mu).grp1(ALL_30ms(mu).grp1(:,size(ALL_30ms(mu).grp1,2))==suj, 1: size(ALL_30ms(mu).grp1,2)-1 );
    Y2= ALL_30ms(mu).grp2(ALL_30ms(mu).grp2(:,size(ALL_30ms(mu).grp2,2))==suj, 1: size(ALL_30ms(mu).grp2,2)-1 );
    Y3= ALL_30ms(mu).grp3(ALL_30ms(mu).grp3(:,size(ALL_30ms(mu).grp3,2))==suj, 1: size(ALL_30ms(mu).grp3,2)-1 );
    Y4= ALL_30ms(mu).grp4(ALL_30ms(mu).grp4(:,size(ALL_30ms(mu).grp4,2))==suj, 1: size(ALL_30ms(mu).grp4,2)-1 );


    A = [zeros(size(Y2,1),1);ones(size(Y4,1),1)];%;2*ones(size(Y3,1),1);3*ones(size(Y4,1),1)];
    subplot(6,2,suj)
    F = spm1d.stats.nonparam.anova1([Y2; Y4],A); % Y1 Y3
    Fi =F.inference(0.05,'iterations', 30);
    Fi.plot()
    %Fi.plot_threshold_label();
    %Fi.plot_p_values();
    
end

%% ANOVA RT

% 9 11
% [ 4, 8 ]
% [2, 3, 5:7, 10]

for mu = [2, 3, 5:7, 10]
    
    Y1= ALL_30ms(mu).grp1;
    Y2= ALL_30ms(mu).grp2;
    Y3= ALL_30ms(mu).grp3;
    Y4= ALL_30ms(mu).grp4;
    Y5= ALL_30ms(mu).grp5;
    Y6= ALL_30ms(mu).grp6;
    Y7= ALL_30ms(mu).grp7;
    Y8= ALL_30ms(mu).grp8;

    sY= [];
    for i = 1:12
        sY(i,1)=size(Y1(Y1(:,size(Y1,2))==i,:),1);
        sY(i,2)=size(Y2(Y2(:,size(Y2,2))==i,:),1);
        sY(i,3)=size(Y3(Y3(:,size(Y3,2))==i,:),1);
        sY(i,4)=size(Y4(Y4(:,size(Y4,2))==i,:),1);
        sY(i,5)=size(Y5(Y5(:,size(Y5,2))==i,:),1);
        sY(i,6)=size(Y6(Y6(:,size(Y6,2))==i,:),1);
        sY(i,7)=size(Y7(Y7(:,size(Y7,2))==i,:),1);
        sY(i,8)=size(Y8(Y8(:,size(Y8,2))==i,:),1);
    end

    ndata=min(min(sY(sY~=0)));
    Yp1 = [];
    Yp2 = [];
    Yp3 = [];
    Yp4 = [];
    Yp5 = [];
    Yp6 = [];
    Yp7 = [];
    Yp8 = [];

    for i = 1:12
        Yp1 = [Yp1; datasample(Y1(Y1(:,size(Y1,2))==i,:),ndata)];
        Yp2 = [Yp2; datasample(Y2(Y2(:,size(Y2,2))==i,:),ndata)];
        Yp3 = [Yp3; datasample(Y3(Y3(:,size(Y3,2))==i,:),ndata)];
        Yp4 = [Yp4; datasample(Y4(Y4(:,size(Y4,2))==i,:),ndata)];
        Yp5 = [Yp5; datasample(Y5(Y5(:,size(Y5,2))==i,:),ndata)];
        Yp6 = [Yp6; datasample(Y6(Y6(:,size(Y6,2))==i,:),ndata)];
        Yp7 = [Yp7; datasample(Y7(Y7(:,size(Y7,2))==i,:),ndata)];
        Yp8 = [Yp8; datasample(Y8(Y8(:,size(Y8,2))==i,:),ndata)];
    end

    Y = [Yp1; Yp2; Yp3; Yp4];% Yp5; Yp6; Yp7; Yp8];
    %%
    h=figure;
    title(ALL_30ms(mu).Names)

    A = [zeros(size(Yp1,1),1); ones(size(Yp2,1),1); 2*ones(size(Yp3,1),1); 3*ones(size(Yp4,1),1)];


    iterations = 200;

    F  = spm1d.stats.nonparam.anova1rm(Y(:,1:size(Y,2)-1), A, Y(:,size(Y,2))); %)
    Fi = F.inference(0.05,'iterations', iterations);
    Fi.plot()
    Fi.plot_threshold_label();
    Fi.plot_p_values();
    
    saveas(h,['C:\Users\p1218107\Documents\Data_Piano\graphiques\',ALL_30ms(mu).Names],'jpeg')
end

