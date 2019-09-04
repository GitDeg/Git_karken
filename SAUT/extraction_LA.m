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
for suj= 2: 12 %2:length(sujet)


repertoireLA =    [pathname,sujet{suj},'\EnregistrementLA']; % Changer pour le nom du répertoire avec tes essais au piano

filesLA = fullfile(repertoireLA, '*.c3d'); 

c3dfilesLA = dir(filesLA); 

%% Création des structures de stockage

CINE_LA = [] ; 
CINE_LA.RawData.Label = [] ;
CINE_LA.RawData.Data = [] ;
CINE_LA.Rate = [] ;
CINE_LA.Filename = [] ;


mark = find(cellfun(@isempty,strfind({c3dfilesLA.name}, 'FraSta')') ==0 );

for m = 1 : length(mark)
    c3dfilesLA_int(m).name = c3dfilesLA(mark(m)).name;
end
%% Ouverture des Data et stockage dans structures

for i= 1 : length(c3dfilesLA_int) % [1 2 12 13]
    
    FileNamesLA = fullfile(repertoireLA, c3dfilesLA_int(i).name);
    c3dDataLA = ezc3dRead(FileNamesLA);
    
    AnalogNames = c3dDataLA.parameters.ANALOG.LABELS; 
    
    count=0;
    
    for j = 1:size(c3dDataLA.data.points,2)
        CINE_LA(i).RawData(j).Label = c3dDataLA.parameters.POINT.LABELS{j,:};  
        CINE_LA(i).RawData(j).Data = permute(c3dDataLA.data.points(:,j,:), [1,3,2]); 
    end
    CINE_LA(i).Rate = c3dDataLA.parameters.POINT.RATE;
    CINE_LA(i).Filename = c3dfilesLA_int(i).name(1,1:end-4);
end

save(['\\10.89.24.15\j\Valentin\SAUT' sujet{suj} '\CINE_LA.mat'], 'CINE_LA')

toc
end