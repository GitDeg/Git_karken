%% Initialisation

clear all;
close all;
clc;


%addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))

pathname='C:\Users\p1218107\Documents\Data_Piano\data_1.5s';

%% Ouverture des datas 


cd(pathname)

%load([pathname,'\All_grp.mat'])
load([pathname,'\Post_ttt_semiauto.mat'])

%%
close all

F=2100/4;


t = 1/F -788/(2*F) : 1/F : 788/(2*F);

for mu = 10 % : 10
    
    sujets = [];
    clear('sY','Y','Y1','Y2','Y3','Y4','Yp1','Yp2','Yp3','Yp4','A','B')
    for grp = 1 : 8
        suff = num2str(grp);
%        assignin('base', sprintf('Y%s',suff), ALL_grp(grp).mu(mu).mean)
%        assignin('base', sprintf('Y%s',suff), [(Post_ttt(grp).suj_mean(1).data+ Post_ttt(grp).suj_mean(2).data)/2 , [1:12]'])
        assignin('base', sprintf('Y%s',suff), [Post_ttt(grp).suj_mean(mu).data , [1:12]'])
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
    
    

    Ypre = [Y1; Y3; Y5; Y7];    % StaTen    = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); ones(size(Y3,1),1); ones(size(Y4,1),1)];
    Yfra = [Y2; Y4; Y6; Y8];    % MsBas     = [zeros(size(Y1,1),1); ones(size(Y2,1),1); zeros(size(Y3,1),1); ones(size(Y4,1),1)];
    
    
    Ysta = [Y1; Y2; Y3; Y4];    % PreFra    = [zeros(size(Y1,1),1); ones(size(Y2,1),1); zeros(size(Y3,1),1); ones(size(Y4,1),1)];
    Yten = [Y5; Y6; Y7; Y8];     MsBas     = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); ones(size(Y3,1),1); ones(size(Y4,1),1)];
    
    
    Yms  = [Y1; Y2; Y5; Y6];    % PreFra    = [zeros(size(Y1,1),1); ones(size(Y2,1),1); zeros(size(Y3,1),1); ones(size(Y4,1),1)];
    Ybas = [Y3; Y4; Y7; Y8];     StaTen    = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); ones(size(Y3,1),1); ones(size(Y4,1),1)];
    
    
    Subj = [suj_set ;suj_set ;suj_set ;suj_set] ;

    iterations = 1500 ; %1000;
    
%     A = [zeros(size(Y1,1),1); ones(size(Y2,1),1); zeros(size(Y3,1),1); ones(size(Y4,1),1); ...
%         zeros(size(Y5,1),1); ones(size(Y6,1),1); zeros(size(Y7,1),1); ones(size(Y8,1),1)]; % pressé (0) / frappé (1)
%     
%     B = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); ones(size(Y3,1),1); ones(size(Y4,1),1); ...
%         zeros(size(Y5,1),1); zeros(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % MS (0) / Bassin (1)
%     
%     C = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); zeros(size(Y3,1),1); zeros(size(Y4,1),1); ...
%         ones(size(Y5,1),1); ones(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % Sta (0) / Tenue (1)
    
%     Subj = [suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set]; %Y(:,size(Y,2));

    figure(mu)    
    suptitle([Post_ttt(8).suj_mean(mu).name(1:end-1) ' n=' num2str(length(suj_set))])
    set(gcf, 'Position', get(0, 'Screensize'))  
%%

    F = spm1d.stats.nonparam.anova1rm(Ysta, MsBas, Subj);
    Fi =F.inference(0.05,'iteration',iterations);

    subplot(2,2,1)
    plot(t, Fi.z,'k')
    hold on 
    line(xlim , [Fi.zstar Fi.zstar] ,'Color',[1 0 0], 'LineStyle','--')
    title('MS / Bassin (Stacato)')
    plot(t, 3*mean(Ysta(logical(MsBas),:)), 'color' ,'c', 'linewidth',2)   %% Bassin
    plot(t, 3*mean(Ysta(~logical(MsBas),:)), 'color' ,'m', 'linewidth',2)  %% MS
    
    i = 1;
    delta = 50;
    while i ~= size(Fi.clusters,2)+1 

        SP = round(Fi.clusters{1, i}.endpoints);
        if SP(1) == 0
            SP(1) = 1;
        end

        while i ~= size(Fi.clusters,2) && abs(SP(2) - round(Fi.clusters{1, i+1}.endpoints(1))) < delta
            SP(2) = round(Fi.clusters{1, i+1}.endpoints(2));

            if SP(2) == round(Fi.clusters{1, size(Fi.clusters,2)}.endpoints(2))
                i = size(Fi.clusters,2);
            else
                i = i+1;
            end
        end
        lim = ylim;
        surf(t(SP(1):SP(2)), [lim(1) : (lim(2) - lim(1))/(SP(2)-SP(1)) : lim(2)], repmat(Fi.z(SP(1):SP(2)), length(Fi.z(SP(1):SP(2))) ,1),'FaceAlpha',0.3, 'EdgeColor','none')
        view(2)
        i = i+1;
    end
    legend({'p value' 'seuil' 'Bassin' 'MS'})
    
%%    
    F = spm1d.stats.nonparam.anova1rm(Yten, MsBas, Subj);
    Fi =F.inference(0.05,'iteration',iterations);

    subplot(2,2,2)
    plot(t, Fi.z,'k')
    hold on 
    line(xlim , [Fi.zstar Fi.zstar] ,'Color',[1 0 0], 'LineStyle','--')
    title('Ms/Bassin (Tenue)')
    plot(t, 3*mean(Yten(logical(MsBas),:)), 'color' ,'c', 'linewidth',2)   %% Bassin
    plot(t, 3*mean(Yten(~logical(MsBas),:)), 'color' ,'m', 'linewidth',2)  %% MS
    
    i = 1;
    delta = 50;
    while i ~= size(Fi.clusters,2)+1 

        SP = round(Fi.clusters{1, i}.endpoints);
        if SP(1) == 0
            SP(1) = 1;
        end

        while i ~= size(Fi.clusters,2) && abs(SP(2) - round(Fi.clusters{1, i+1}.endpoints(1))) < delta
            SP(2) = round(Fi.clusters{1, i+1}.endpoints(2));

            if SP(2) == round(Fi.clusters{1, size(Fi.clusters,2)}.endpoints(2))
                i = size(Fi.clusters,2);
            else
                i = i+1;
            end
        end
        lim = ylim;
        surf(t(SP(1):SP(2)), [lim(1) : (lim(2) - lim(1))/(SP(2)-SP(1)) : lim(2)], repmat(Fi.z(SP(1):SP(2)), length(Fi.z(SP(1):SP(2))) ,1),'FaceAlpha',0.3, 'EdgeColor','none')
        view(2)
        i = i+1;
    end
    legend({'p value' 'seuil' 'Bassin' 'MS'})
    
%%    
    F = spm1d.stats.nonparam.anova1rm(Yms, StaTen, Subj);
    Fi =F.inference(0.05,'iteration',iterations);

    subplot(2,2,3)
    plot(t, Fi.z,'k')
    hold on 
    line(xlim , [Fi.zstar Fi.zstar] ,'Color',[1 0 0], 'LineStyle','--')
    title('Sta / Ten  (MS)')
    plot(t, 3*mean(Ypre(logical(StaTen),:)), 'color' ,'b', 'linewidth',2)   %% Ten
    plot(t, 3*mean(Ypre(~logical(StaTen),:)), 'color' ,'y', 'linewidth',2)  %% Staccato
    
    i = 1;
    delta = 50;
    while i ~= size(Fi.clusters,2)+1 

        SP = round(Fi.clusters{1, i}.endpoints);
        if SP(1) == 0
            SP(1) = 1;
        end

        while i ~= size(Fi.clusters,2) && abs(SP(2) - round(Fi.clusters{1, i+1}.endpoints(1))) < delta
            SP(2) = round(Fi.clusters{1, i+1}.endpoints(2));

            if SP(2) == round(Fi.clusters{1, size(Fi.clusters,2)}.endpoints(2))
                i = size(Fi.clusters,2);
            else
                i = i+1;
            end
        end
        lim = ylim;
        surf(t(SP(1):SP(2)), [lim(1) : (lim(2) - lim(1))/(SP(2)-SP(1)) : lim(2)], repmat(Fi.z(SP(1):SP(2)), length(Fi.z(SP(1):SP(2))) ,1),'FaceAlpha',0.3, 'EdgeColor','none')
        view(2)
        i = i+1;
    end
    legend({'p value' 'seuil' 'Tenue' 'Staccato'})
    
%%    
    F = spm1d.stats.nonparam.anova1rm(Ybas, StaTen, Subj);
    Fi =F.inference(0.05,'iteration',iterations);

    subplot(2,2,4)
    plot(t, Fi.z,'k')
    hold on 
    line(xlim , [Fi.zstar Fi.zstar] ,'Color',[1 0 0], 'LineStyle','--')
    title('Sta / Ten (Bassin)')
    plot(t, 3*mean(Yfra(logical(StaTen),:)), 'color' ,'b', 'linewidth',2)   %% Ten
    plot(t, 3*mean(Yfra(~logical(StaTen),:)), 'color' ,'y', 'linewidth',2)  %% Staccato
    
    i = 1;
    delta = 50;
    while i ~= size(Fi.clusters,2)+1 

        SP = round(Fi.clusters{1, i}.endpoints);
        if SP(1) == 0
            SP(1) = 1;
        end

        while i ~= size(Fi.clusters,2) && abs(SP(2) - round(Fi.clusters{1, i+1}.endpoints(1))) < delta
            SP(2) = round(Fi.clusters{1, i+1}.endpoints(2));

            if SP(2) == round(Fi.clusters{1, size(Fi.clusters,2)}.endpoints(2))
                i = size(Fi.clusters,2);
            else
                i = i+1;
            end
        end
        lim = ylim;
        surf(t(SP(1):SP(2)), [lim(1) : (lim(2) - lim(1))/(SP(2)-SP(1)) : lim(2)], repmat(Fi.z(SP(1):SP(2)), length(Fi.z(SP(1):SP(2))) ,1),'FaceAlpha',0.3, 'EdgeColor','none')
        view(2)
        i = i+1;
    end
    legend({'p value' 'seuil' 'Tenue' 'Staccato'})
  %%  
    saveas(figure(1), '\\10.89.24.15\e\Bureau\Valentin\graphiques\chronologie', 'jpeg')
    
end

