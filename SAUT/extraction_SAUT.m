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
%% Cr�ation des chemins d'acces

pathname = '\\10.89.24.15\e\Projet_piano\Bassin';


%%
tic
for suj= 5%2:length(sujet)

repertoireSAUT =    [pathname,sujet{suj},'\EnregistrementSAUTS']; % Changer pour le nom du r�pertoire avec tes essais au piano

filesSAUT = fullfile(repertoireSAUT, '*.c3d'); 

c3dfilesSAUT = dir(filesSAUT); 

%% Cr�ation des structures de stockage
EMG_SAUT=[];
EMG_SAUT.RawData=[];
EMG_SAUT.Labels={};
EMG_SAUT.Rate=[];
EMG_SAUT.Filename=[];
EMG_SAUT.reps = [];

CINE_SAUT = [] ; 
CINE_SAUT.RawData.Label = [] ;
CINE_SAUT.RawData.Data = [] ;
CINE_SAUT.Rate = [] ;
CINE_SAUT.Filename = [] ;

mark1 = find(cellfun(@isempty,strfind({c3dfilesSAUT.name}, 'Ant60Aig')') ==0 );
mark2 = find(cellfun(@isempty,strfind({c3dfilesSAUT.name}, 'AntBas60Aig')') ==0 );
mark3 = find(cellfun(@isempty,strfind({c3dfilesSAUT.name}, 'Dem60Aig')') ==0 );
mark4 = find(cellfun(@isempty,strfind({c3dfilesSAUT.name}, 'DemBas60Aig')') ==0 );

c3dfilesSAUT_int(1).name = c3dfilesSAUT(mark1).name;
c3dfilesSAUT_int(2).name = c3dfilesSAUT(mark2).name;
c3dfilesSAUT_int(3).name = c3dfilesSAUT(mark3).name;
c3dfilesSAUT_int(4).name = c3dfilesSAUT(mark4).name;

%% Ouverture des Data EMG et stockage dans structures

for i=1:length(c3dfilesSAUT_int)
    
    FileNamesSAUT = fullfile(repertoireSAUT, c3dfilesSAUT_int(i).name);
    c3dDataSAUT = ezc3dRead(FileNamesSAUT);
    
    EMG_Vicon_Name = {'Triceps.EMG1','Biceps.EMG2',...
                      'DeltAnt.EMG3','DeltMed.EMG4',...
                      'TrapSup.EMG5','GrandDent.EMG6',...
                      'GranDent.EMG6','GrandPec.IM EMG7',...
                      'Extenseurs.EMG8','Extenseur.EMG8','Sensor 8.EMG8',...
                      'Flechisseurs1.EMG9','Flechisseurs2.EMG10','FlechisseursA.EMG9','FlechisseursB.EMG10',...
                      'Fléchisseurs1.EMG9','Fléchisseurs2.EMG10','Fléchisseur1.EMG9',...
                      'Sensor 9.EMG9','Sensor 10.EMG10'};
                  
    AnalogNames = c3dDataSAUT.parameters.ANALOG.LABELS; 
    
    count=0;
    for h=1:size(EMG_Vicon_Name,2)
        [im, wh]= ismember(EMG_Vicon_Name{h}, AnalogNames);
        if im==1
            count=count+1; 
            EMG_SAUT(i).RawData(count,:) = c3dDataSAUT.data.analogs(:,wh)';  %Matrice EMG, fichiers LA
            EMG_SAUT(i).Labels    = [EMG_SAUT(i).Labels EMG_Vicon_Name(h)] ; %Nom des entr�es EMG
        end
    end
    EMG_SAUT(i).Rate = c3dDataSAUT.parameters.ANALOG.RATE; %fr�quence d'�chantillonage.
    EMG_SAUT(i).Filename = c3dfilesSAUT_int(i).name(1,1:end-4); %Nom de l'essai
    
    for j = 1:size(c3dDataSAUT.data.points,2)
        CINE_SAUT(i).RawData(j).Label = c3dDataSAUT.parameters.POINT.LABELS{j,:};  
        CINE_SAUT(i).RawData(j).Data = permute(c3dDataSAUT.data.points(:,j,:), [1,3,2]); 
    end
    CINE_SAUT(i).Rate = c3dDataSAUT.parameters.POINT.RATE;
    CINE_SAUT(i).Filename = c3dfilesSAUT_int(i).name(1,1:end-4);
end
end
%% Sauvegarde

% save(['C:\Users\p1218107\Documents\SAUT' sujet{suj} '\EMG_SAUT'],'EMG_SAUT');
save(['C:\Users\p1218107\Documents\SAUT' sujet{suj} '\CINE_SAUT'],'CINE_SAUT');
% 
% save(['\\10.89.24.15\j\Valentin\SAUT' sujet{suj} '\EMG_SAUT'],'EMG_SAUT');
save(['\\10.89.24.15\j\Valentin\SAUT' sujet{suj} '\CINE_SAUT'],'CINE_SAUT');
toc
% end