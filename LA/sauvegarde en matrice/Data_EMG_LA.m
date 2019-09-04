%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%Chargement des fichier c3d des EMG pour un sujet / exportation du fichier%                                                                        %
%                              Protocole LA                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation

clear all;
close all;
clc;

addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))

sujet = {'\001','\002','\003','\004','\005','\006','\007','\008','\009','\010','\011','\012'}; % Changer le nom du sujet
%% Création des chemins d'acces

pathname = 'C:\Users\Deg\Desktop\Thèse\DATA_Piano';

%%
tic
for j=1:length(sujet)


repertoireLA =    [pathname,sujet{j},'\LA']; % Changer pour le nom du répertoire avec tes essais au piano

filesLA = fullfile(repertoireLA, '*.c3d'); 

c3dfilesLA = dir(filesLA); 

%% Création des structures de stockage
EMG_LA=[];
EMG_LA.RawData=[];
EMG_LA.Labels={};
EMG_LA.Rate=[];
EMG_LA.Filename=[];
EMG_LA.reps = [];


%% Ouverture des Data EMG et stockage dans structures

for i=1:length(c3dfilesLA)
    
    FileNamesLA = fullfile(repertoireLA, c3dfilesLA(i).name);
    c3dDataLA = ezc3dRead(FileNamesLA);
    
    EMG_Vicon_Name = {'Triceps.EMG1','Biceps.EMG2',...
                      'DeltAnt.EMG3','DeltMed.EMG4',...
                      'TrapSup.EMG5','GrandDent.EMG6',...
                      'GranDent.EMG6','GrandPec.IM EMG7',...
                      'Extenseurs.EMG8','Extenseur.EMG8','Sensor 8.EMG8',...
                      'Flechisseurs1.EMG9','Flechisseurs2.EMG10','FlechisseursA.EMG9','FlechisseursB.EMG10',...
                      'FlÃ©chisseurs1.EMG9','FlÃ©chisseurs2.EMG10','FlÃ©chisseur1.EMG9',...
                      'Sensor 9.EMG9','Sensor 10.EMG10'};
                  
    AnalogNames = c3dDataLA.parameters.ANALOG.LABELS; 
    
    count=0;
    for h=1:size(EMG_Vicon_Name,2)
        [im, wh]= ismember(EMG_Vicon_Name{h}, AnalogNames);
        if im==1
            count=count+1; 
            EMG_LA(i).RawData(count,:) = c3dDataLA.data.analogs(:,wh)'; %Matrice EMG, fichiers LA
            EMG_LA(i).Labels    = [EMG_LA(i).Labels EMG_Vicon_Name(h)] ; %Nom des entrées EMG
        end
    end
    EMG_LA(i).Rate = c3dDataLA.parameters.ANALOG.RATE; %fréquence d'échantillonage.
    EMG_LA(i).Filename = c3dfilesLA(i).name(1,1:end-4); %Nom de l'essai

end
toc

%% Sauvegarde

save([repertoireLA,'\all_EMG_',sujet{j}(2:end)],'EMG_LA');
end