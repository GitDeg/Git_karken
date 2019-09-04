%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Chargement des fichier c3d pour un sujet / exportation du fichier force %                                                                        %
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
Fo_LA=[];
Fo_LA.RawData=[];
Fo_LA.Labels={};
Fo_LA.Rate=[];
Fo_LA.Filename=[];
Fo_LA.reps = []; % ?

%% Ouverture des Data EMG et stockage dans structures

for i=1:length(c3dfilesLA) 
    
    FileNamesLA = fullfile(repertoireLA, c3dfilesLA(i).name);
    c3dDataLA = ezc3dRead(FileNamesLA);
    
    Force_Vicon_Name = {'CapteurForceGauche.EMG11', 'CapteurForceDroit.EMG12',...
                        'CapteurDeForceGauche.EMG11', 'CapteurDeForceDroit.EMG12'};
                  
    AnalogNames = c3dDataLA.parameters.ANALOG.LABELS; 
    
    count=0;
    for h=1:size(Force_Vicon_Name,2)
        [im, wh]= ismember(Force_Vicon_Name{h}, AnalogNames);
        if im==1
            count=count+1; 
            Fo_LA(i).RawData(count,:) = c3dDataLA.data.analogs(:,wh)'-mean(c3dDataLA.data.analogs(1:150,wh)'); %Matrice EMG, fichiers LA
            Fo_LA(i).Labels    = [Fo_LA(i).Labels Force_Vicon_Name(h)] ; %Nom des entrées EMG
        end
    end
    Fo_LA(i).Rate = c3dDataLA.parameters.ANALOG.RATE; %fréquence d'échantillonage.
    Fo_LA(i).Filename = c3dfilesLA(i).name(1,1:end-4); %Nom de l'essai

end
toc

%% Sauvegarde

save([repertoireLA,'\all_Fo_',sujet{j}(2:end)],'Fo_LA');
end