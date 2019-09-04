%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                           Découpage cycles                        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc 

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
%%
% pathname = '\\10.89.24.15\j\Valentin\SAUT';
pathname = 'C:\Users\p1218107\Documents\SAUT';
sujet = {'\001','\002','\003','\004','\005','\006','\007','\008','\009','\010','\011','\012'};
load('C:\Users\p1218107\Documents\SAUT\allREPS.mat')

GrpNames = {' AntiHoraire ', ...        %% 60 aigue toujours
            ' AntiHoraire Bassin ', ...
            ' Demicercle ', ...
            ' Demicercle Bassin '...
            ' Frappé Staccato '...
            ' Frappé Staccato Bassin'};
%% Sélection des Datas à traiter

for suj = 1 : 12

    load([pathname sujet{suj} '\CINE_SAUT.mat'])
    load([pathname sujet{suj} '\EMG_SAUT.mat'])
    load([pathname sujet{suj} '\CINE_LA.mat'])
    load([pathname sujet{suj} '\EMG_LA.mat'])
    Data = [];
    
    for grpi = 1 : length(CINE_SAUT)
        
        
            
        MuscleNames = {'Triceps','Biceps','DeltAnt','DeltMed',...   % Uniformisation des noms des muscles
                   'TrapSup','GrandDent','GrandPec','Extenseurs','Fléchisseurs'};

        Data(grpi).Filename = GrpNames{grpi};
        
        DATA = EMG_SAUT(grpi).RawData;
        
        %% préparation du filtrage
        Fs = EMG_SAUT(grpi).Rate; %Féquence d'échantillonnage
        Ts = 1/Fs;
        
        [b,a] = butter(2,[5 500]/(Fs/2),'bandpass'); %Filtre passe bande, fenetre recommandée
        wind_length = round(30*10^(-3)/Ts); % fenetre de 30ms pour filtrage rms
        
        %% filtrage         
        
        tp = DATA  - mean(DATA,2) ;% Soustraction baseline
        
        tp_filt = filtfilt(b,a,tp')'; % Application filtre passe bande
                 
        [yupper,ylower] = envelope(tp_filt',wind_length,'rms');
        Data(grpi).EMG.EMG_env = yupper';
        
        %% Cine
        
        mark1 = find(cellfun(@isempty,strfind({CINE_SAUT(grpi).RawData.Label}, 'Meta1')') ==0 );
        for i = 1 : length(mark1)
            Data(grpi).CINE(i).name = CINE_SAUT(grpi).RawData(mark1(i)).Label;
            Data(grpi).CINE(i).data = CINE_SAUT(grpi).RawData(mark1(i)).Data;  
            Data(grpi).CINE(i).data(Data(grpi).CINE(i).data == 0) = NaN;
        end
        
        mark2 = find(cellfun(@isempty,strfind({CINE_SAUT(grpi).RawData.Label}, 'meta5')') ==0 );
        for i = length(mark1)+1 : (length(mark1) + length(mark2))
            Data(grpi).CINE(i).name = CINE_SAUT(grpi).RawData(mark2(i-length(mark1))).Label;
            Data(grpi).CINE(i).data = CINE_SAUT(grpi).RawData(mark2(i-length(mark1))).Data;  
            Data(grpi).CINE(i).data(Data(grpi).CINE(i).data == 0) = NaN;
        end
        
        
%%
        
        for cycle = 2 : length(EMG_SAUT(grpi).reps)-2
            
            interv_EMG = round(EMG_SAUT(grpi).reps(cycle)*Fs) : round(EMG_SAUT(grpi).reps(cycle + 1)*Fs);
            interv_CINE = round(EMG_SAUT(grpi).reps(cycle)*150) : round(EMG_SAUT(grpi).reps(cycle + 1)*150);
            
            t = 50;
            p = 150;

            for mu = 1 : 8
                Data(grpi).EMG.mu(mu).name = MuscleNames{mu};
                Data(grpi).EMG.mu(mu).cycle(cycle,:) = resample(Data(grpi).EMG.EMG_env(mu,interv_EMG),1000, length(Data(grpi).EMG.EMG_env(mu,interv_EMG))); 
            end
            Data(grpi).EMG.mu(9).name = MuscleNames{9};
            Data(grpi).EMG.mu(9).cycle(cycle,:) = resample((Data(grpi).EMG.EMG_env(9,interv_EMG) + Data(grpi).EMG.EMG_env(10,interv_EMG))/2 ,1000, length(Data(grpi).EMG.EMG_env(9,interv_EMG))); 
            
            for mark = 1 : length(Data(grpi).CINE)
                 q_CINE = length(Data(grpi).CINE(mark).data(:,interv_CINE(1)-t : interv_CINE(end)+t));
                 A = resample(Data(grpi).CINE(mark).data(:,interv_CINE(1)-t : interv_CINE(end)+t)', p, q_CINE)';
                 Data(grpi).CINE(mark).cycle(cycle*3-2 :cycle*3 ,:) = A(:, length(A)/2 - 55 : length(A)/2 + 55);
            end
            
        end
    end    
   
    %% Partie LA 
    reps_i = find(cellfun(@isempty,strfind({allREPS(suj).sujet.name}, 'FraSta')') ==0 );
    
    for grpla = 5 : 6
        Data(grpla).Filename = EMG_LA(grpla-4).name;
        
        for mu = 1 : 9
            Data(grpla).EMG.mu(mu).name = MuscleNames{mu};
            Data(grpla).EMG.mu(mu).cycle = resample(EMG_LA(grpla-4).mu(mu).cycle'  ,1000, length(EMG_LA(grpla-4).mu(mu).cycle'))'; 
        end
        
        mark1 = find(cellfun(@isempty,strfind({CINE_LA( (grpla-4)*2 -1 ).RawData.Label}, 'meta3')') ==0 );
        mark2 = find(cellfun(@isempty,strfind({CINE_LA( (grpla-4)*2 ).RawData.Label}, 'meta3')') ==0 );
        
        for i = 1 : length(mark1)
            Data(grpla).CINE(i).name = CINE_LA( (grpla-4)*2 -1 ).RawData(mark1(i)).Label;
            Data(grpla).CINE(i).data1 = CINE_LA( (grpla-4)*2 -1 ).RawData(mark1(i)).Data; 
            Data(grpla).CINE(i).data1(Data(grpla).CINE(i).data1 == 0) = NaN;
            Data(grpla).CINE(i).reps1 = allREPS(suj).sujet(reps_i( (grpla-4)*2 -1 )).true_REPS - allREPS(suj).sujet(reps_i( (grpla-4)*2 -1 )).first_frame/150;
            
            Data(grpla).CINE(i).data2 = CINE_LA( (grpla-4)*2 ).RawData(mark2(i)).Data;
            Data(grpla).CINE(i).data2(Data(grpla).CINE(i).data2 == 0) = NaN;
            Data(grpla).CINE(i).reps2 = allREPS(suj).sujet(reps_i( (grpla-4)*2 )).true_REPS - allREPS(suj).sujet(reps_i( (grpla-4)*2 )).first_frame/150;
            
            for cyc = 2 : length(Data(grpla).CINE(i).reps1)-5 
                interv1 = round( (Data(grpla).CINE(i).reps1(cyc) - 0.75)*150) : round((Data(grpla).CINE(i).reps1(cyc) + 0.75)*150);
                Data(grpla).CINE(i).cycle((cyc-1)*3-2 : (cyc-1)*3,:) = Data(grpla).CINE(i).data1(:,interv1);
            end
            
            ncycles = size(Data(grpla).CINE(i).cycle,1);
            
            for cyc = 2 : length(Data(grpla).CINE(i).reps2)-5 
                interv2 = round( (Data(grpla).CINE(i).reps2(cyc) - 0.75)*150) : round((Data(grpla).CINE(i).reps2(cyc) + 0.75)*150);
                Data(grpla).CINE(i).cycle(ncycles + (cyc-1)*3-2 : ncycles + (cyc-1)*3,:) = Data(grpla).CINE(i).data2(:,interv2);
            end
        end
    end
    
%% Choix 10 cycles + normalisation

    for mu = 1 : 9
        for grpi = 1 : length(Data)
            RMSE = [];
            for i = 1 : size(Data(grpi).EMG.mu(mu).cycle,1)
                RMSE(i) = sqrt( mean( (Data(grpi).EMG.mu(mu).cycle(i,:) - mean(Data(grpi).EMG.mu(mu).cycle)).^2 )); 
            end 
            [s,i]=sort(RMSE);
            
            Data(grpi).EMG.mu(mu).cycle_select = Data(grpi).EMG.mu(mu).cycle(i(1:10), :); 
        end

%% Normalisation
        
        if suj == 1 && mu == 6
            n_grp = [1 2 4 5 6];
            
        elseif  suj == 2 && mu == 7
            n_grp = [3 4 5 6];
            
        elseif suj == 6 && mu == 8
            n_grp = [1 2 5 6];
            
        else
            n_grp = 1 : length(Data);
        end
        
        all_acti = [];
        for grpi = n_grp
            all_acti = [all_acti, Data(grpi).EMG.mu(mu).cycle_select];
        end
        
        m_acti = mean(mean(all_acti));
        
        for grpi = 1 : 6
            Data(grpi).EMG.mu(mu).norm = Data(grpi).EMG.mu(mu).cycle_select / m_acti; 
        end 
%%
    end
    
    for grpi = 1 : length(Data)
        for mark = 1 : length(Data(grpi).CINE)
            M_Cine = [mean(Data(grpi).CINE(mark).cycle(1:3:end,:));...
                    mean(Data(grpi).CINE(mark).cycle(2:3:end,:));...
                    mean(Data(grpi).CINE(mark).cycle(3:3:end,:))];
                
            RMSE_cine = []; 
            
            for i = 1 : size(Data(grpi).CINE(mark).cycle,1)/3

                RMSE_X(i) = sqrt( mean( (Data(grpi).CINE(mark).cycle(i*3-2 ,:) - M_Cine(1,:) ).^2 ));
                RMSE_Y(i) = sqrt( mean( (Data(grpi).CINE(mark).cycle(i*3-1 ,:) - M_Cine(2,:) ).^2 ));
                RMSE_Z(i) = sqrt( mean( (Data(grpi).CINE(mark).cycle(i*3 ,:) - M_Cine(3,:) ).^2 ));

                RMSE_cine(i) = sqrt( RMSE_X(i).^2 + RMSE_Y(i).^2 + RMSE_Z(i).^2);
            end

            [s,i]=sort(RMSE_cine);
            i_CINE = [i(1:10)*3-2;    i(1:10)*3-1;    i(1:10)*3];

            Data(grpi).CINE(mark).cycle_select = Data(grpi).CINE(mark).cycle(i_CINE(:), :);

        end
    end
     
    save([pathname sujet{suj} '\CYCLES.mat'], 'Data')
end

%%
for suj = 1 : 12
    load([pathname sujet{suj} '\CYCLES.mat'])
    
    for grpi = 1 : 6
       All_data(grpi).Filename = GrpNames(grpi);
        
       for mu =  1 : 9
            All_data(grpi).EMG(suj).mu(mu).name  = Data(grpi).EMG.mu(mu).name;
            All_data(grpi).EMG(suj).mu(mu).data  = Data(grpi).EMG.mu(mu).norm;
         
            All_data(grpi).EMG(1).mean_EMG(mu).name = Data(grpi).EMG.mu(mu).name;
            All_data(grpi).EMG(1).mean_EMG(mu).data(suj,:) = mean(Data(grpi).EMG.mu(mu).norm);
       end
       
       for mark = 1 : length(Data(grpi).CINE)
            All_data(grpi).CINE(suj).mark(mark).name = Data(grpi).CINE(mark).name;  
            All_data(grpi).CINE(suj).mark(mark).data = Data(grpi).CINE(mark).cycle_select;
            
            All_data(grpi).CINE(suj).mark(mark).mean_data(1,:) = nanmean(Data(grpi).CINE(mark).cycle_select((1:10)*3-2,:));
            All_data(grpi).CINE(suj).mark(mark).mean_data(2,:) = nanmean(Data(grpi).CINE(mark).cycle_select((1:10)*3-1,:));
            All_data(grpi).CINE(suj).mark(mark).mean_data(3,:) = nanmean(Data(grpi).CINE(mark).cycle_select((1:10)*3,:));
       end
    end
     
end

All_data(1).EMG(2).mu(7).data = [];
All_data(1).EMG(1).mean_EMG(7).data(2,:) = NaN;   

All_data(2).EMG(2).mu(7).data = [];
All_data(2).EMG(1).mean_EMG(7).data(2,:) = NaN;  

All_data(3).EMG(1).mu(6).data = [];
All_data(3).EMG(1).mean_EMG(6).data(1,:) = NaN;  

All_data(3).EMG(6).mu(8).data = [];
All_data(3).EMG(1).mean_EMG(8).data(6,:) = NaN;  

All_data(4).EMG(6).mu(8).data = [];
All_data(4).EMG(1).mean_EMG(8).data(6,:) = NaN;  

save([pathname '\All_data.mat'], 'All_data')

%%

grp = 6

for suj = 1 : 12
    figure(suj)
    subplot(3,1,1)
    hold on
    plot(All_data(grp).CINE(suj).mark.data((1:10)*3 -2,:)')
    subplot(3,1,2)
    hold on
    plot(All_data(grp).CINE(suj).mark.data((1:10)*3-1,:)')
    subplot(3,1,3)
    hold on
    plot(All_data(grp).CINE(suj).mark.data((1:10)*3,:)')
end