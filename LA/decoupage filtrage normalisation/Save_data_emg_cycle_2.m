%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Sauvegarde des cycles de chaque groupe et chaque sujet        %
%                              Protocole LA                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation

clear all;
close all;
clc;

%addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB\spm1dmatlab-master'))
sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};

pathname='C:\Users\p1218107\Documents\Data_Piano\data_1.5s'; 


%% Sélection des Datas à traiter

for n= 4 %1:length(sujet)
        
%         load([pathname,'\',sujet{n},'_EMG_GRP.mat'])
        
        PathName_EMG=[pathname(1:end-10),'\',sujet{n},'\'];
        
        load([PathName_EMG,'All_Grp1.mat'])
        load([PathName_EMG,'All_Grp2.mat'])
        load([PathName_EMG,'All_Grp3.mat'])
        load([PathName_EMG,'All_Grp4.mat'])             
        load([PathName_EMG,'All_Grp5.mat'])
        load([PathName_EMG,'All_Grp6.mat'])
        load([PathName_EMG,'All_Grp7.mat'])
        load([PathName_EMG,'All_Grp8.mat'])
        
        Grps=[];            % création d'une matrice contenant l'ensemble des datas
        Grps{1}=Grp1.EMG;
        Grps{2}=Grp2.EMG;
        Grps{3}=Grp3.EMG;
        Grps{4}=Grp4.EMG;
        Grps{5}=Grp5.EMG;
        Grps{6}=Grp6.EMG;
        Grps{7}=Grp7.EMG;
        Grps{8}=Grp8.EMG;        


        
        GRP = [] ;          % création d'une matrice à sauvegarder
        GRP.Name = [] ;
        GRP.EMG = [];       % EMG non aligné
        GRP.EMG.cycle = [] ;
        GRP.EMG.cycle_filt = [] ;
        GRP.EMG.cycle_env = [] ;
        GRP.EMG.mean = [] ;
        GRP.EMG.labels= [] ; 
        
        GRP.EMGal = [];     % EMG aligné
        GRP.EMGal.cycle = [] ;
        GRP.EMGal.mean = [] ;
        GRP.EMGal.labels= [] ;
        
        
    %% Groupe 1: Membre supérieur pressé staccato -> Grp1
    %% Groupe 2: Membre supérieur frappé staccato -> Grp2
    %% Groupe 3: Bassin pressé staccato -> Grp3
    %% Groupe 4: Bassin frappé staccato -> Grp4
    %% Groupe 5: Membre supérieur pressé tenue ->Grp5
    %% Groupe 6: Membre supérieur frappé tenue ->Grp6
    %% Groupe 7: Bassin pressé tenue ->Grp7
    %% Groupe 8: Bassin frappé tenue ->Grp8
    
    GrpNames={'Groupe 1: Membre supérieur pressé staccato',...
              'Groupe 2: Membre supérieur frappé staccato',...
              'Groupe 3: Bassin pressé staccato',...
              'Groupe 4: Bassin frappé staccato',...
              'Groupe 5: Membre supérieur pressé tenue',...
              'Groupe 6: Membre supérieur frappé tenue',...
              'Groupe 7: Bassin pressé tenue'...
              'Groupe 8: Bassin frappé tenue'};      
          
    
    for j = 5:8
        Grps{1,j}(1).Labels={'Triceps','Biceps','DeltAnt','DeltMed',...   % Uniformisation des noms des muscles
                            'TrapSup','GrandDent','GrandPec','Extenseurs','Fléchisseurs1','Fléchisseurs2'};
        Fs = Grps{1,j}(1).Rate; %Féquence d'échantillonnage
        Ts = 1/Fs;
        
        nbframes = round(Fs + Fs/2); % taille de la fenêtre d'analyse
        P=2100; Q=Fs;
        
        [b,a]=butter(2,[20 300]/(Fs/2),'bandpass'); %Filtre passe bande


        for k=1:size(Grps{1,j}(1).RawData,1)
            
            ki=[9 10 8 2 1 3 4 7 5 6]; % rangement des datas en fonction de la distance entre le muscle et la touche
            
            emg = [];
            emg_filt = [];
            emg_env = [];
            emg_resamp = [];
            tp=[];
            tp_filt=[];
            

            
            for p=1:2
                
                tp = Grps{1,j}(p).RawData(ki(k),:)-mean(Grps{1,j}(p).RawData(ki(k),:));% Soustraction baseline
                tp_filt = filtfilt(b,a,tp); % Application filtre passe bande
                
                wl = round(30*10^(-3)/Ts); 
                
                [yupper,ylower] = envelope(tp,wl,'rms');
                
                REPS = [] ; 
                REPS = Grps{1,j}(p).reps(2:end-1);

                for m = 1 : length(REPS)
                    
                    tdeb = REPS(m)-floor(nbframes/2);
                    tfin = REPS(m)+floor(nbframes/2);
                    
                    emg(m+length(Grps{1,j}(1).reps(2:end-1))*(p-1),:)=tp(tdeb+1:tfin);
                    emg_filt(m+length(Grps{1,j}(1).reps(2:end-1))*(p-1),:)=tp_filt(tdeb+1:tfin);
                    emg_env(m+length(Grps{1,j}(1).reps(2:end-1))*(p-1),:)=yupper(tdeb+1:tfin);
                    emg_resamp(m+length(Grps{1,j}(1).reps(2:end-1))*(p-1),:)=resample(yupper(tdeb+1:tfin),P,Q);
                   
                    
                end
            end
          
            Mean_emg = mean(emg,1);
            Mean_emg_filt = mean(emg_filt,1);
            Mean_emg_env = mean(emg_env,1);
           
            emg_ali = [];
            emg_ali_filt = [];
            emg_ali_env = [];
            emg_ali_resamp = [];
            tp=[];
            tp_filt=[];
            
            D = [];
            D_filt = [];
            D_env = [];
            
            for p=1:2
                
                tp = Grps{1,j}(p).RawData(ki(k),:)-mean(Grps{1,j}(p).RawData(ki(k),:)); % Soustraction baseline
               
                tp_filt = filtfilt(b,a,tp); % Application filtre passe bande
                
                wl = round(30*10^(-3)/Ts); 
                
                [yupper,ylower] = envelope(tp,wl,'rms');
                REPS = [];
                REPS = Grps{1,j}(p).reps(2:end-1);
                
                for m=1:length(REPS)

                    tdeb = REPS(m)-floor(nbframes/2);
                    tfin = REPS(m)+floor(nbframes/2);

                    delta=100;
                    
                    D(m) = finddelay(Mean_emg,tp(tdeb+delta:tfin-delta));
                    D_filt(m) = finddelay(Mean_emg_filt,tp_filt(tdeb+delta:tfin-delta));
                    D_env(m) = finddelay(Mean_emg_env,yupper(tdeb+delta:tfin-delta));

                    
                    emg_ali(m+length(Grps{1,j}(1).reps(2:end-1))*(p-1),:)=tp(tdeb+1+delta+D(m):tfin+delta+D(m));
                    emg_ali_filt(m+length(Grps{1,j}(1).reps(2:end-1))*(p-1),:)=tp_filt(tdeb+1+delta+D_filt(m):tfin+delta+D_filt(m));
                    emg_ali_env(m+length(Grps{1,j}(1).reps(2:end-1))*(p-1),:)=yupper(tdeb+1+delta+D_env(m):tfin+delta+D_env(m));
                    emg_ali_resamp(m+length(Grps{1,j}(1).reps(2:end-1))*(p-1),:)=resample(yupper(tdeb+1+delta+D_env(m):tfin+delta+D_env(m)),P,Q);
                
                end
            end
% 
            Mean_emg_ali_env=mean(emg_ali_env,1);

            GRP(j).Name = GrpNames{j};
            
            GRP(j).EMG(k).cycle = emg;
            GRP(j).EMG(k).cycle_filt = emg_filt;
            GRP(j).EMG(k).cycle_env = emg_env;
            GRP(j).EMG(k).cycle_resamp = emg_resamp;
            GRP(j).EMG(k).mean_env = Mean_emg_env;
            GRP(j).EMG(k).labels = Grps{1,j}(1).Labels{ki(k)};
            GRP(j).Rate=2100;

            GRP(j).EMGal(k).cycle = emg_ali;
            GRP(j).EMGal(k).cycle_ali_filt = emg_ali_filt;
            GRP(j).EMGal(k).cycle_ali_env = emg_ali_env;
            GRP(j).EMGal(k).cycle_ali_resamp = emg_ali_resamp;
            GRP(j).EMGal(k).mean_env = Mean_emg_ali_env;
            GRP(j).EMGal(k).labels = Grps{1,j}(1).Labels{ki(k)};
            
            
        end
    end
    %save([pathname,'\',sujet{n},'_EMG_GRP.mat'],'GRP')
end


