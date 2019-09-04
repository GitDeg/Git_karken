
clear all;
close all;
clc;

%addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB\MoCapWorkshopMontreal2019\Workshop_Material\MocapToolbox_v1.5.1\mocaptoolbox'))
sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};

pathname='C:\Users\p1218107\Documents\Data_Piano\'; 
pathname_Cine = '\\10.89.24.15\e\Projet_Reconstructions\Piano_Reconstructions\script_matlab\';
pathname_c3d = '\\10.89.24.15\e\Projet_piano\Bassin\';
%%

for suj = 9
    
    load([pathname, sujet{suj}, '\NewREPS.mat'])
    load([pathname, 'allREPS.mat'])
    
    clear 'Data' 'GRP'
     
    for grp = 1 : 16
        Data(grp).name = allREPS(suj).sujet(grp).name;
    end
        
    n = {allREPS(suj).sujet.name};
    
    for grp = 1 : 8
        for i = 1 : 2
            place = find(strcmp(n,GRPs(grp).EMG(i).Filename));
            
            Data(place).EMG = GRPs(grp).EMG(i).RawData ;
            Data(place).Rate = GRPs(grp).EMG(i).Rate ;
            
        end
    end
    
    
    
    for grp = 1 : 16 
        
        t = 1/Data(grp).Rate: 1/Data(grp).Rate: length(Data(grp).EMG)/Data(grp).Rate;
        
        figure(grp)
%         figure(suj)
%         subplot(4, 4, grp)
        plot(t, Data(grp).EMG(10,:), 'k')
        hold on 
    %    line([Data(9).REPS;Data(9).REPS]*Data(9).Rate, repmat(ylim, length(Data(9).REPS),1)'  )
        line([allREPS(suj).sujet(grp).true_REPS; allREPS(suj).sujet(grp).true_REPS ] , repmat(ylim, length(allREPS(suj).sujet(grp).true_REPS),1)','Color', 'r')
        line([allREPS(suj).sujet(grp).reps_fo; allREPS(suj).sujet(grp).reps_fo ]/allREPS(suj).sujet(grp).rate.force , repmat(ylim, length(allREPS(suj).sujet(grp).reps_fo),1)', 'Color', 'b')
        xlim([0 t(length(t))])
        
    end
    
end


%%

suj = 12; % sujet 9 : [1 2 3 9 10 11]

load([pathname_Cine 'allDataPiano_' sujet{suj}])
load([pathname, 'allREPS.mat'])

for grp = setdiff(1:16, [8 13 15])

    c3d_cine = readc3d([pathname_c3d sujet{suj} '\EnregistrementLA\' allREPS(suj).sujet(grp).name '.c3d']);
    FirstFrame_cine = c3d_cine.Header.FirstVideoFrame ; 

    c3d_fo = readc3d(['\\10.89.24.15\e\Projet_piano\Bassin\Data_Temporaire\' sujet{suj} '\' allREPS(suj).sujet(grp).name '.c3d']);
    FirstFrame_fo = c3d_fo.Header.FirstVideoFrame ; 

    frame_dif = FirstFrame_cine - FirstFrame_fo; 

    a = allData(suj,grp).Tobs(:,68,:);
    data = a(:,:);
    sdata = sort(data); 
    Cine = -(data(3,:));
    Cine([1:200, length(Cine)-300:length(Cine)]) = NaN;

    t = 1/Data(grp).Rate: 1/Data(grp).Rate: length(Data(grp).EMG)/Data(grp).Rate;
%   t_cine = 1/150 : 1/150 : length(Cine)/150;
%   frame_dif = (t(length(t)) - t_cine(length(t_cine)))*150;

    allREPS(suj).sujet(grp).first_frame = frame_dif;

    Cine = [NaN(1, frame_dif), Cine];
    allREPS(suj).sujet(grp).reps_cine = allREPS(suj).sujet(grp).reps_cine + frame_dif;

%     Cine = [NaN(1, round(frame_dif)), Cine]; Cine=Cine(abs(round(frame_dif)):end);
%     allREPS(suj).sujet(grp).reps_cine = allREPS(suj).sujet(grp).reps_cine + round(frame_dif);

    Cine_rep = allREPS(suj).sujet(grp).reps_cine(1) /allREPS(suj).sujet(grp).rate.cine;
    Fo = allREPS(suj).sujet(grp).reps_fo /allREPS(suj).sujet(grp).rate.force;
    Son = allREPS(suj).sujet(grp).reps_son /allREPS(suj).sujet(grp).rate.son;
    Son(Son == 0) = [];

    [a b] = min(abs(Fo - Cine_rep));
    Son = Son - (Son(b) - Cine_rep);
    allREPS(suj).sujet(grp).true_REPS = Son;

    t_cine = 1/150 : 1/150 : length(Cine)/150;
    figure

    plot(t, 250*Data(grp).EMG(10,:), 'k')
    hold on 
    plot(t_cine, Cine, 'g')

    line([allREPS(suj).sujet(grp).true_REPS; allREPS(suj).sujet(grp).true_REPS ] , repmat(ylim, length(allREPS(suj).sujet(grp).true_REPS),1)','Color', 'r')
    line([Fo; Fo ], repmat(ylim, length(allREPS(suj).sujet(grp).reps_fo),1)', 'Color', 'b')
    line([allREPS(suj).sujet(grp).reps_cine; allREPS(suj).sujet(grp).reps_cine ]/allREPS(suj).sujet(grp).rate.cine , repmat(ylim, length(allREPS(suj).sujet(grp).reps_cine),1)', 'Color', 'm')
    xlim([0 t(length(t))])
end


%%
save([pathname 'allREPS.mat'],'allREPS')

%%
load([pathname, 'allREPS.mat'])

for grp = 1 : 16
    for suj = [1:8 10:12]
        
        
        n = length(allREPS(suj).sujet(grp).reps_fo);
        m = length(allREPS(suj).sujet(grp).true_REPS);
        
        if n > 20
            n = 20;
        end
        
        Fo = allREPS(suj).sujet(grp).reps_fo /allREPS(suj).sujet(grp).rate.force;
        True_reps = allREPS(suj).sujet(grp).true_REPS;

        clear 'reps_i' 'min_dif' 'dif'

        for i = 1 : n
            reps_i = ones(1,m)*Fo(i);
            dif = diff([reps_i; True_reps]);
            min_dif(i) = dif( find( abs(dif) == min( abs( dif))));
        end

        Dif_Fo_Cine(grp,suj) = nanmean(min_dif);
    end
end

Dif_Fo_Cine(:,9) = [];

%%
boxplot(Dif_Fo_Cine', [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8])
hold on
line( xlim, [mean(mean(Dif_Fo_Cine)) , mean(mean(Dif_Fo_Cine))], 'Color', 'm')

%%

for grp = 1 : 16
    allREPS(9).sujet(grp).true_REPS = (allREPS(9).sujet(grp).reps_fo + round(mean(mean(Dif_Fo_Cine))*allREPS(9).sujet(grp).rate.force))/allREPS(9).sujet(grp).rate.force  ;
end

