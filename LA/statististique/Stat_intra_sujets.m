%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     Analyse statistique intrasujet                      %
%                              Protocole LA                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation

clear all;
close all;
clc;

%addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB\spm1dmatlab-master'))

sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};

pathname = 'C:\Users\p1218107\Documents\Data_Piano';
pathname_graph = 'C:\Users\p1218107\Documents\Data_Piano\graphiques\statistique par sujet\ANOVA2\';
F=2100;


t = 1/F-525/F: 1/F : 524/(F);
%% Ouverture des datas 


for n = [3:8 , 10:12]    
    close all
    pathname_emg = ['C:\Users\p1218107\Documents\Data_Piano\',sujet{n},'\EMG_grp_30ms.mat'];
    
    cd(pathname)
    
    load(pathname_emg)

%save([pathname,'\All_data_muscle_30ms.mat'],'ALL_30ms')
%%
    for mu = 1:10
%% ANOVA2 + MEAN_SD cloud

        Y1= GRP(1).EMGal(mu).cycl_norm;
        Y2= GRP(2).EMGal(mu).cycl_norm;
        Y3= GRP(3).EMGal(mu).cycl_norm;
        Y4= GRP(4).EMGal(mu).cycl_norm;

        sY(1)=size(Y1,1);
        sY(2)=size(Y2,1);
        sY(3)=size(Y3,1);
        sY(4)=size(Y4,1);

        ndata=min(sY);
        Yp1 = [];
        Yp2 = [];
        Yp3 = [];
        Yp4 = [];
        
        Yp1 = [Yp1; datasample(Y1,ndata)];
        Yp2 = [Yp2; datasample(Y2,ndata)];
        Yp3 = [Yp3; datasample(Y3,ndata)];
        Yp4 = [Yp4; datasample(Y4,ndata)];

        Y = [Yp1; Yp2; Yp3; Yp4];
        
        A = [zeros(size(Yp1,1),1);ones(size(Yp2,1),1);zeros(size(Yp3,1),1);ones(size(Yp4,1),1)];
        B = [zeros(size(Yp1,1),1);zeros(size(Yp2,1),1);ones(size(Yp3,1),1);ones(size(Yp4,1),1)];
        
        F = spm1d.stats.nonparam.anova2(Y,A,B);
        Fi =F.inference(0.01,'iterations', 500);
        
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
        suptitle(['sujet ', sujet{n},' ', GRP(1).EMGal(mu).labels]) 
        
        %legend(GRP(1).Name,'', GRP(2).Name,'', GRP(3).Name,'', GRP(4).Name)
        legend('MS pressé','', 'MS frappé','', 'Bassin pressé','', 'Bassin frappé')

        xlim([0 1049])
        
        saveas(figure(mu),[pathname_graph, 'sujet ', sujet{n},' ', GRP(1).EMGal(mu).labels],'jpeg')
    end
end

