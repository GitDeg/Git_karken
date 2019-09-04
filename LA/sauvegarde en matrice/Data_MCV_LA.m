%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%Chargement des fichier c3d des MCV pour un sujet / exportation du fichier%                                                                        %
%                              Protocole LA                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisation

clear all;
close all;
clc;

addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))

sujet = {'\001','\002','\003','\004','\005','\006','\007','\008','\009','\010','\011','\012'}; % Changer le nom du sujet


%%
%% Chargement de tous les fichiers

pathname='C:\Users\Deg\Desktop\Thèse\DATA_Piano';

for n=2%1:length(sujet)
repertoireMVC =     [pathname,sujet{n},'/MVC']; % Changer pour le nom du répertoire avec tes MVC
filesMVC = fullfile(repertoireMVC, '*.c3d'); 
c3dfilesMVC = dir(filesMVC); 


%% Ouverture de tous les fichiers MVC
% Création d'une structure contenant seulement les données EMG de chaque fichier 

FileName = ['C:\Users\Deg\Desktop\Thèse\DATA_Piano',sujet{n},'\MVC','\EMG_MVC.mat'];


% if exist(FileName) && TagMVC
%     load(FileName)
% else

    MVC=[];
    EMG_MVC.RawData=[];
    EMG_MVC.Labels=[];
    EMG_MVC.Rate=[];
    EMG_MVC.Filename=[];
%%

    for i=[1:18 20]%length(c3dfilesMVC)
        FileNamesMVC = fullfile(repertoireMVC, c3dfilesMVC(i).name);
        c3dDataMVC = ezc3dRead(FileNamesMVC);

        % Selection du data du fichier c3d correspondant à l'EMG
        EMG_Vicon_Name = {'Triceps.EMG1','Biceps.EMG2','DeltAnt.EMG3','DeltMed.EMG4','TrapSup.EMG5','GrandDent.EMG6','GranDent.EMG6','GrandPec.IM EMG7','Extenseurs.EMG8','Extenseur.EMG8','Sensor 8.EMG8','Flechisseurs1.EMG9','Flechisseurs2.EMG10','FlechisseursA.EMG9','FlechisseursB.EMG10','FlÃ©chisseurs1.EMG9','FlÃ©chisseurs2.EMG10','FlÃ©chisseur1.EMG9','Sensor 9.EMG9','Sensor 10.EMG10'};


        for h=1:20
            for p=1:length(c3dDataMVC.parameters.ANALOG.LABELS)
                if strcmp(c3dDataMVC.parameters.ANALOG.LABELS{p},EMG_Vicon_Name{h})
                    
                    break
                end    
            end
            % Création de la structure EMG_MVC(RawData,Labels,Rate)
            c3dDataMVC.data.analogs_dim1 = size(c3dDataMVC.data.analogs, 1);
            MVC(i).RawData(1:c3dDataMVC.data.analogs_dim1,h) = c3dDataMVC.data.analogs(:,p); %Matrice des EMG des fichiers MVC
            MVC(i).Labels(h,1) = c3dDataMVC.parameters.ANALOG.LABELS(p,1); %nom des EMG selon l'ordre EMG_Vicon_Name.
            EMG_MVC(i).RawData = MVC(i).RawData';% Matrice MVC renversée : lignes (10 EMG dans l'ordre EMG_Vicon_Name), colonnes (ensemble de frames de l'essai). 
            EMG_MVC(i).Labels = MVC(i).Labels;

        end
        IdxMVC_0 = find(EMG_MVC(i).RawData(:,1)==0); %trouver les lignes en trop = 0
        EMG_MVC(i).RawData(IdxMVC_0,:)=[]; %supression des lignes en trop.
        EMG_MVC(i).Labels(IdxMVC_0)=[]; %supression des lignes en trop 
        EMG_MVC(i).Rate = c3dDataMVC.parameters.ANALOG.RATE; %fréquence d'échantillonage.
        EMG_MVC(i).Filename = c3dfilesMVC(i).name(1,1:end-4); %Nom de l'essai

    end
    
    save(FileName, 'EMG_MVC');
    
    %clear('EMG_MVC','MVC','IdxMVC_0')
end


