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

pathname='C:\Users\p1218107\Documents\Data_Piano\'; 


%% Sélection des Datas à traiter

for suj = 10 %1:length(sujet) 
        
    load([pathname, sujet{suj}, '\NewREPS.mat'])
    load([pathname, 'allREPS.mat'])
     
%% préparation des datas
    
    clear 'Data' 'GRP'
     
    for grp = 1 : 16
%        Data(grp).name = allREPS(suj).sujet(grp).name;
%        Data(grp).REPS = allREPS(suj).sujet(grp).true_REPS;
        Data(grp).REPS_corr = allREPS(suj).sujet(grp).REPS_corr;
    end
        
    n = {allREPS(suj).sujet.name};
    
    for grp = 1 : 8
        for i = 1 : 2
            place = find(strcmp(n,GRPs(grp).EMG(i).Filename));
            
            Data(place).EMG = GRPs(grp).EMG(i).RawData ;
            Data(place).Rate = GRPs(grp).EMG(i).Rate ;
            
        end
    end

    GrpNames = {' Bassin frappé staccato ', ...
                ' Bassin frappé tenue ', ...
                ' Bassin pressé staccato ', ...
                ' Bassin pressé tenue ', ...
                ' Membre supérieur frappé staccato ', ...
                ' Membre supérieur frappé tenue ', ...
                ' Membre supérieur pressé staccato ', ...
                ' Membre supérieur pressé tenue '};
            
    MuscleNames = {'Triceps','Biceps','DeltAnt','DeltMed',...   % Uniformisation des noms des muscles
                   'TrapSup','GrandDent','GrandPec','Extenseurs','Fléchisseurs1','Fléchisseurs2' };
          
%% 
    for grp = 1 : 16
%% préparation du filtrage
        Fs = Data(grp).Rate; %Féquence d'échantillonnage
        Ts = 1/Fs;
        
         nbframes = round(Fs + Fs/2); % taille de la fenêtre d'analyse
         P=2100; Q=Fs; % variables pour le resample
        
        
        [b,a] = butter(2,[5 500]/(Fs/2),'bandpass'); %Filtre passe bande, fenetre recommandée
        wind_length = round(30*10^(-3)/Ts); % fenetre de 30ms pour filtrage rms
        
%% filtrage         
        
        tp = Data(grp).EMG  - mean(Data(grp).EMG,2) ;% Soustraction baseline
        
        Data(grp).EMG_filt = filtfilt(b,a,tp')'; % Application filtre passe bande
                 
        [yupper,ylower] = envelope(Data(grp).EMG_filt',wind_length,'rms');
        Data(grp).EMG_env = yupper';
         
%% Découpage 

        REPS = floor(Data(grp).REPS * Fs);
%        REPS = floor(Data(grp).REPS_corr * Fs);
        
        tdeb = REPS - floor(nbframes/2);
        tfin = REPS + floor(nbframes/2);
        
        for mu = 1 : 10
            for cycle = 2 : length(REPS)-1
                
%                Data(grp).cycle(mu).muscle = MuscleNames{mu};
                Data(grp).cycle(mu).EMG(cycle,:) = Data(grp).EMG_env(mu, tdeb(cycle) : tfin(cycle) );
                Data(grp).cycle(mu).EMG_resamp(cycle,:) = resample(Data(grp).EMG_env(mu, tdeb(cycle)+1 : tfin(cycle) ), P, Q);
            
            end
        end
        
    end

    save([pathname,'\data_newREPS\',sujet{suj},'.mat'],'Data')
    
    for grp = 1 : 8
        grpbis = grp*2-1;
        
        GRP(grp).name = GrpNames{grp};
        GRP(grp).rate = 2100;
        
        for mu = 1 : 8
            
            GRP(grp).cycle(mu).muscle = MuscleNames{mu};
        
            GRP(grp).cycle(mu).data = [ Data(grpbis).cycle(mu).EMG_resamp; Data(grpbis + 1).cycle(mu).EMG_resamp ];
            
            GRP(grp).cycle(mu).data(GRP(grp).cycle(mu).data(:,1)==0,:) = [] ;
        end
        
        GRP(grp).cycle(9).muscle = 'Flechisseurs';
        GRP(grp).cycle(9).data = [(Data(grpbis).cycle(9).EMG_resamp + Data(grpbis).cycle(10).EMG_resamp)/2 ;...
                                    (Data(grpbis+1).cycle(9).EMG_resamp + Data(grpbis+1).cycle(10).EMG_resamp)/2];
                                
        GRP(grp).cycle(9).data(GRP(grp).cycle(9).data(:,1)==0,:) = [] ;                        
    end
    
    save([pathname,'\data_newREPS\',sujet{suj},'_data.mat'],'GRP')
end


