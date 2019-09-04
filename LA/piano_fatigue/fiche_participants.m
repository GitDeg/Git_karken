clear all;
close all;
clc;

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
cd('\\10.89.24.15\j\Piano_Fatigue')

%% Fiche participants

path_parti = 'Fiches_participants\travaildirig\';


for suj = 1 : 50
    
    
    if suj <10
        S_suj = sprintf('S10%s', num2str(suj));
    else
        S_suj = sprintf('S1%s', num2str(suj));
    end
    
    [num,txt,raw] = xlsread([path_parti S_suj '.xlsx']);
    
    
    Info_participants(suj).Etiquette     = num(1,1);
    Info_participants(suj).Age           = num(1,2);
    Info_participants(suj).Sexe          = num(1,3);  %% 0 = homme et 1 = femme
    Info_participants(suj).Poids         = num(1,4);
    Info_participants(suj).Lateralite    = num(1,5);  %% 0 = droitier et 1 = gaucher
    Info_participants(suj).nb_annees     = num(1,6);
    Info_participants(suj).nb_heures     = num(1,7);
    Info_participants(suj).MVC           = [num(1,8) num(1,9)];
    Info_participants(suj).Doigtes       = num(1,10);
    Info_participants(suj).Crt_Palm      = num(1,11);
    Info_participants(suj).Hauteur       = num(1,12); %% 0 = absent et 1 = present
    Info_participants(suj).Av_bras       = num(1,13);
    Info_participants(suj).Ordre         = raw(2,14);
    Info_participants(suj).Hanon         = num(:,15);
    Info_participants(suj).Liszt         = num(:,16);
    Info_participants(suj).Temps         = num(:,17);
end

%%
