%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                         D�coupage cycles                        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc 

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
%%
% pathname = '\\10.89.24.15\j\Valentin\SAUT';
pathname = 'C:\Users\p1218107\Documents\Data_Piano\SAUT';
sujet = {'\001','\002','\003','\004','\005','\006','\007','\008','\009','\010','\011','\012'};

GrpNames = {' AntiHoraire ', ...        %% 60 aigue toujours
            ' AntiHoraire Bassin ', ...
            ' Demicercle ', ...
            ' Demicercle Bassin '...
            ' Frapp� Staccato '...
            ' Frapp� Staccato Bassin'};
%% S�lection des Datas � traiter

for suj = 5 %11: 12

    load([pathname sujet{suj} '\CINE_SAUT.mat'])
    load([pathname sujet{suj} '\CINE_LA.mat'])
    Piano = [];
    
    for grpi = 1 : 2 %length(CINE_SAUT)
        %% Cine
        Piano(grpi).GrpNames = GrpNames{grpi};
        
        mark1 = find(cellfun(@isempty,strfind({CINE_SAUT(grpi).RawData.Label}, 'GraDown')') ==0 );
        mark2 = find(cellfun(@isempty,strfind({CINE_SAUT(grpi).RawData.Label}, 'AigDown')') ==0 );
        mark3 = find(cellfun(@isempty,strfind({CINE_SAUT(grpi).RawData.Label}, 'AigAnt')') ==0 );
        
        Piano(grpi).CINE(1).name = CINE_SAUT(grpi).RawData(mark1).Label;
        Piano(grpi).CINE(1).data = CINE_SAUT(grpi).RawData(mark1).Data;  
        Piano(grpi).CINE(1).data(Piano(grpi).CINE(1).data == 0) = NaN;

        Piano(grpi).CINE(2).name = CINE_SAUT(grpi).RawData(mark2).Label;
        Piano(grpi).CINE(2).data = CINE_SAUT(grpi).RawData(mark2).Data;  
        Piano(grpi).CINE(2).data(Piano(grpi).CINE(2).data == 0) = NaN;
        
        Piano(grpi).CINE(3).name = CINE_SAUT(grpi).RawData(mark3).Label;
        Piano(grpi).CINE(3).data = CINE_SAUT(grpi).RawData(mark3).Data;  
        Piano(grpi).CINE(3).data(Piano(grpi).CINE(3).data == 0) = NaN;

    end    
   
    %% Partie LA 
    for grpla = 5 : 6
        
        Piano(grpla).GrpNames = GrpNames{grpla};
        
        mark1 = find(cellfun(@isempty,strfind({CINE_LA((grpla-4)*2 -1).RawData.Label}, 'GraDown')') ==0 );
        mark2 = find(cellfun(@isempty,strfind({CINE_LA((grpla-4)*2 -1).RawData.Label}, 'AigDown')') ==0 );
        mark3 = find(cellfun(@isempty,strfind({CINE_LA((grpla-4)*2 -1).RawData.Label}, 'AigAnt')') ==0 );
        
        mark12 = find(cellfun(@isempty,strfind({CINE_LA((grpla-4)*2).RawData.Label}, 'GraDown')') ==0 );
        mark22 = find(cellfun(@isempty,strfind({CINE_LA((grpla-4)*2).RawData.Label}, 'AigDown')') ==0 );
        mark32 = find(cellfun(@isempty,strfind({CINE_LA((grpla-4)*2).RawData.Label}, 'AigAnt')') ==0 );
        
        Piano(grpla).CINE(1).name = CINE_LA((grpla-4)*2 -1).RawData(mark1).Label;
        Piano(grpla).CINE(1).data = [CINE_LA((grpla-4)*2 -1).RawData(mark1).Data, CINE_LA((grpla-4)*2).RawData(mark12).Data];  
        Piano(grpla).CINE(1).data(Piano(grpla).CINE(1).data == 0) = NaN;

        Piano(grpla).CINE(2).name = CINE_LA((grpla-4)*2 -1).RawData(mark2).Label;
        Piano(grpla).CINE(2).data = [CINE_LA((grpla-4)*2 -1).RawData(mark2).Data, CINE_LA((grpla-4)*2).RawData(mark22).Data];  
        Piano(grpla).CINE(2).data(Piano(grpla).CINE(2).data == 0) = NaN;
        
        Piano(grpla).CINE(3).name = CINE_LA((grpla-4)*2 -1).RawData(mark3).Label;
        Piano(grpla).CINE(3).data = [CINE_LA((grpla-4)*2 -1).RawData(mark3).Data, CINE_LA((grpla-4)*2).RawData(mark32).Data];  
        Piano(grpla).CINE(3).data(Piano(grpla).CINE(3).data == 0) = NaN;
        
    end
    
    save([pathname sujet{suj} '\PIANO.mat'], 'Piano')
end

%%
clear all
close all
clc 

%% Ubuntu 

% addpath(genpath('/media/valentin/060C6AFC0C6AE5E1/Users/p1218107/Documents/MATLAB'))
% 
% pathname = '/media/valentin/060C6AFC0C6AE5E1/Users/p1218107/Documents/Data_Piano/SAUT';
% sujet = {'/001','/002','/003','/004','/005','/006','/007','/008','/009','/010','/011','/012'};
% 
% load([pathname '/All_data1s.mat'])

%% Windows
% pathname = '\\10.89.24.15\j\Valentin\SAUT';
addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
addpath(genpath('\\10.89.24.15\e\Librairies\S2M_Lib\Benjamin_dev'))

pathname = 'C:\Users\p1218107\Documents\Data_Piano\SAUT';

sujet = {'\001','\002','\003','\004','\005','\006','\007','\008','\009','\010','\011','\012'};

load([pathname '\All_data1s_xy.mat'])

%%
X = (1:3:30);
Y = (1:3:30) + 1;
Z = (1:3:30) + 2;

for suj = 1 : 12
    load([pathname sujet{suj} '\PIANO.mat'])
    Gra_down = [];
    Aig_down = [];
    Aig_ant  = [];
    
    if suj == 3
        grp_i = [1 2 6];
    elseif suj == 8
        grp_i = [1 5 6]; 
    else
        grp_i = [1 2 5 6];
    end
    
    for grp = grp_i
        Gra_down = [Gra_down , Piano(grp).CINE(1).data];
        Aig_down = [Aig_down , Piano(grp).CINE(2).data];
        Aig_ant  = [Aig_ant  , Piano(grp).CINE(3).data];
    end
    Gra_down = nanmean(Gra_down')';
    Aig_down = nanmean(Aig_down')';
    Aig_ant  = nanmean(Aig_ant')';
    
   
    X_axis = (Aig_down - Gra_down)/ norm(Aig_down - Gra_down) ;
    v = (Aig_ant - Aig_down)/ norm(Aig_ant - Aig_down);
    Z_axis = cross(X_axis,v)/ norm(cross(X_axis,v));
    Y_axis = cross(Z_axis, X_axis);
    P_G_P = [X_axis Y_axis Z_axis Gra_down; 0 0 0 1];
  
    P_P_T = [[0.999997943568855,-0.00183443587290552,-0.000864698264544730,846.441959231975];...
             [0.00173439311223167,0.552601595008462,0.833443740797563,20.2973130215416];...
             [-0.00105106545597917,-0.833443526604612,0.552603640254274,-43.9616852345589];...
             [0,0,0,1]];
    
    P_G_T = P_G_P*P_P_T; 
    
    A(:,suj) = invR(P_G_P)*[Aig_down;1]; 
    G(:,suj) = invR(P_G_T)*[Gra_down;1];
    d_piano(suj) = norm(A(:,suj));
    
    %%
    for grp = [1 2 5 6]
        for mark = 1 : length(All_data(grp).CINE(suj).mark)
            All_data(grp).CINE(suj).mark(mark).data_transp = nan(size(All_data(grp).CINE(suj).mark(mark).data));
            
            for i = 1 : 10
                cycl = invR(P_G_T)*[All_data(grp).CINE(suj).mark(mark).data(i*3-2:i*3,:); ones(1,length(All_data(grp).CINE(suj).mark(mark).data))];
                All_data(grp).CINE(suj).mark(mark).data_transp(i*3-2:i*3,:) = cycl(1:3,:);
            end
        end
    end
end
