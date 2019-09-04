%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                           Découpage cycles                        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc 

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
%%
pathname = '\\10.89.24.15\j\Valentin\SAUT';
sujet = {'\001','\002','\003','\004','\005','\006','\007','\008','\009','\010','\011','\012'};
load('\\10.89.24.15\j\Valentin\allREPS.mat')
%%

for suj = 1 : 12

    load([pathname sujet{suj} '\CINE_LA.mat'])
    %%
    ii = 0;
    for i = 1 : length(CINE_LA)
        if strcmp({CINE_LA(i).Filename(4:9)}, 'FraSta')
            ii = ii + 1;
            Data_LA(ii).name = CINE_LA(i).Filename;
            mark = find(cellfun(@isempty,strfind({CINE_LA(i).RawData.Label}, 'meta2')') ==0 );
%             if isempty(mark)
%                 mark = find(strcmp({CINE_LA(i).RawData.Label}, ['meta2_post'])==1);
%             end
            Data_LA(ii).CINE.name = 'Meta_2_post';
            Data_LA(ii).CINE.data = CINE_LA(i).RawData(mark).Data;  
            Data_LA(ii).CINE.data(Data_LA(ii).CINE(1).data == 0) = NaN;
            
        end
    end
    
    %%
     ii = 0;
     for i = 1 : length(allREPS(suj).sujet)
        if strcmp({allREPS(suj).sujet(i).name(4:9)}, 'FraSta')
            ii = ii + 1;
            Data_LA(ii).CINE.REPS = allREPS(suj).sujet(i).true_REPS;
        end
     end  
     %%
     for i = 1 : length(Data_LA)  
         for cyc = 2 : length(Data_LA(i).CINE.REPS)-5        
             interv = round((Data_LA(i).CINE.REPS(cyc)-0.75)*150) : round((Data_LA(i).CINE.REPS(cyc) + 0.75)*150);
             Data_LA(i).CINE.cycle((cyc-1)*3-2 : (cyc-1)*3,:) = Data_LA(i).CINE.data(:,interv);
         end
     end
    
     save([pathname sujet{suj} '\CYCLE_LA.mat'], 'Data_LA')
end