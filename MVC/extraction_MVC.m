%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%Chargement des fichier c3d des EMG pour un sujet / exportation du fichier%                                                                        %
%                            Protocole SAUT                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation

clear all;
close all;
clc;

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
addpath(genpath('C:\Users\Deg\Documents\MATLAB'))

sujet = {'\001','\002','\003','\004','\005','\006','\007','\008','\009','\010','\011','\012'}; % Changer le nom du sujet
%% Création des chemins d'acces

pathname = '\\10.89.24.15\e\Projet_piano\Bassin';


%%
tic
for suj= 1:length(sujet)

repertoireMVC =    [pathname,sujet{suj},'\MVC']; % Changer pour le nom du répertoire avec tes essais au piano

filesMVC = fullfile(repertoireMVC, '*.c3d'); 

c3dfilesMVC = dir(filesMVC); 

%% Création des structures de stockage
MVC=[];
MVC.RawData=[];
MVC.Labels={};
MVC.Rate=[];
MVC.Filename=[];


%% Ouverture des Data EMG et stockage dans structures

for i=1%:length(c3dfilesMVC)
    
    FileNamesMVC = fullfile(repertoireMVC, c3dfilesMVC(i).name);
    c3dDataMVC = ezc3dRead(FileNamesMVC);
    
    EMG_Vicon_Name = {'Triceps.EMG1','Biceps.EMG2',...
                      'DeltAnt.EMG3','DeltMed.EMG4',...
                      'TrapSup.EMG5','GrandDent.EMG6',...
                      'GranDent.EMG6','GrandPec.IM EMG7',...
                      'Extenseurs.EMG8','Extenseur.EMG8','Sensor 8.EMG8',...
                      'Flechisseurs1.EMG9','Flechisseurs2.EMG10','FlechisseursA.EMG9','FlechisseursB.EMG10',...
                      'FlÃ©chisseurs1.EMG9','FlÃ©chisseurs2.EMG10','FlÃ©chisseur1.EMG9',...
                      'Sensor 9.EMG9','Sensor 10.EMG10'};
                  
    AnalogNames = c3dDataMVC.parameters.ANALOG.LABELS; 
    
    count=0;
    for h=1:size(EMG_Vicon_Name,2)
        [im, wh]= ismember(EMG_Vicon_Name{h}, AnalogNames);
        if im==1
            count=count+1
            MVC(i).RawData(count,:) = c3dDataMVC.data.analogs(:,wh)'; %Matrice EMG, fichiers LA
            MVC(i).Labels    = [MVC(i).Labels EMG_Vicon_Name(h)] ; %Nom des entrées EMG
        end
    end
    MVC(i).Rate = c3dDataMVC.parameters.ANALOG.RATE; %fréquence d'échantillonage.
    MVC(i).Filename = c3dfilesMVC(i).name(1,1:end-4); %Nom de l'essai
    
end

%% Sauvegarde

% save(['C:\Users\p1218107\Documents\SAUT' sujet{suj} '\EMG_SAUT'],'EMG_SAUT');
% save(['C:\Users\p1218107\Documents\SAUT' sujet{suj} '\CINE_SAUT'],'CINE_SAUT');
% 
save(['\\10.89.24.15\j\Valentin\SAUT\' sujet{suj} '\MVC'],'MVC');
toc
end