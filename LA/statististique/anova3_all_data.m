%% Initialisation

clear all;
close all;
clc;


%addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))

pathname='C:\Users\p1218107\Documents\Data_Piano\data_newREPS';

%% Ouverture des datas 


cd(pathname)

load([pathname,'\All_data.mat'])

%%
close all

F=2100/5;


t = 1/F -630/(2*F) : 1/F : 630/(2*F);

for mu = 1 :9
    
    suj_set = [];
    clear('sY','Y','Y1','Y2','Y3','Y4','Yp1','Yp2','Yp3','Yp4','A','B')
    for grp = 1 : 8
        suff = num2str(grp);
%        assignin('base', sprintf('Y%s',suff), ALL_grp(grp).mu(mu).data)
        assignin('base', sprintf('Y%s',suff), [All_data(grp).suj(1).mu(mu).data_selec;...
                                               All_data(grp).suj(2).mu(mu).data_selec;...
                                               All_data(grp).suj(3).mu(mu).data_selec;...
                                               All_data(grp).suj(4).mu(mu).data_selec;...
                                               All_data(grp).suj(5).mu(mu).data_selec;...
                                               All_data(grp).suj(6).mu(mu).data_selec;...
                                               All_data(grp).suj(7).mu(mu).data_selec;...
                                               All_data(grp).suj(8).mu(mu).data_selec;...
                                               All_data(grp).suj(9).mu(mu).data_selec;...
                                               All_data(grp).suj(10).mu(mu).data_selec;...
                                               All_data(grp).suj(11).mu(mu).data_selec;...
                                               All_data(grp).suj(12).mu(mu).data_selec])
    end
    
    suj_set = mintersect(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8);
    
    Y1 = Y1(ismember(Y1(:,end),suj_set),:);
    Y2 = Y2(ismember(Y2(:,end),suj_set),:);
    Y3 = Y3(ismember(Y3(:,end),suj_set),:);
    Y4 = Y4(ismember(Y4(:,end),suj_set),:);
    Y5 = Y5(ismember(Y5(:,end),suj_set),:);
    Y6 = Y6(ismember(Y6(:,end),suj_set),:);
    Y7 = Y7(ismember(Y7(:,end),suj_set),:);
    Y8 = Y8(ismember(Y8(:,end),suj_set),:);
    
    sY= nan(12,8);
    for i = suj_set'
        sY(i,1)=size(Y1(Y1(:,size(Y1,2))==i,:),1);
        sY(i,2)=size(Y2(Y2(:,size(Y2,2))==i,:),1);
        sY(i,3)=size(Y3(Y3(:,size(Y3,2))==i,:),1);
        sY(i,4)=size(Y4(Y4(:,size(Y4,2))==i,:),1);
        sY(i,5)=size(Y5(Y5(:,size(Y5,2))==i,:),1);
        sY(i,6)=size(Y6(Y6(:,size(Y6,2))==i,:),1);
        sY(i,7)=size(Y7(Y7(:,size(Y7,2))==i,:),1);
        sY(i,8)=size(Y8(Y8(:,size(Y8,2))==i,:),1);
    end

    ndata=nanmin(nanmin(sY));
    Yp1 = [];
    Yp2 = [];
    Yp3 = [];
    Yp4 = [];
    Yp5 = [];
    Yp6 = [];
    Yp7 = [];
    Yp8 = [];

    for i = suj_set'
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

    iterations = 500 ; %1000    
    A = [zeros(size(Yp1,1),1); zeros(size(Yp2,1),1); ones(size(Yp3,1),1); ones(size(Yp4,1),1); ...
        zeros(size(Yp5,1),1); zeros(size(Yp6,1),1); ones(size(Yp7,1),1); ones(size(Yp8,1),1)]; % frappé (0) / pressé (1)
    
    B = [zeros(size(Yp1,1),1); zeros(size(Yp2,1),1); zeros(size(Yp3,1),1); zeros(size(Yp4,1),1); ...
        ones(size(Yp5,1),1); ones(size(Yp6,1),1); ones(size(Yp7,1),1); ones(size(Yp8,1),1)]; % Bassin (0) / MS (1)
    
    C = [zeros(size(Yp1,1),1); ones(size(Yp2,1),1); zeros(size(Yp3,1),1); ones(size(Yp4,1),1); ...
        zeros(size(Yp5,1),1); ones(size(Yp6,1),1); zeros(size(Yp7,1),1); ones(size(Yp8,1),1)]; % Sta (0) / Tenue (1)
    
    
    Subj = Y(:,size(Y,2));
    
    F = spm1d.stats.nonparam.anova3rm(Y(:,1:end-1), A, B, C, Subj);
    Fi =F.inference(0.05,'iteration',iterations);
    
    figure(mu)
    set(gcf, 'Position', get(0, 'Screensize'))
    suptitle([All_data(1).suj_mean(mu).name ' n=' num2str(length(suj_set))])
    subplot(3,3,8:9)
    title('Means All grps')
    hold on
    plot(t, mean(Yp1(:,1:end-1)), 'r')
    plot(t, mean(Yp2(:,1:end-1)), 'g')
    plot(t, mean(Yp3(:,1:end-1)), 'b')
    plot(t, mean(Yp4(:,1:end-1)), 'k')
    plot(t, mean(Yp5(:,1:end-1)), 'y')
    plot(t, mean(Yp6(:,1:end-1)), 'm')
    plot(t, mean(Yp7(:,1:end-1)), 'c')
    plot(t, mean(Yp8(:,1:end-1)), 'r')
    
    
    Titles = {'1. Pressé/Frappé' '2. Bassin / Membre Supérieur' '3. Sta/Ten' 'Interact A/B' 'Interact A/C' 'Interact B/C' 'Interact A/B/C'};
    
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
            if SP(2) == 0
                SP(2) = 1;
            end
            subplot(3,3,8:9)
            hold on
        	line(t([SP(1),SP(1)]),ylim,'Color',[1 0 0])
            line(t([SP(2),SP(2)]),ylim,'Color',[1 0 0])
        end
    end

%    saveas(figure(mu), ['C:\Users\p1218107\Documents\Data_Piano\graphiques\ANOVA3\all_data_nonpara_rm\' ALL_grp(1).mu(mu).name], 'jpeg')
    
end

