%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                          Recherche des reps                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc 

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
%%
% pathname = '\\10.89.24.15\j\Valentin\SAUT';
pathname = 'C:\Users\p1218107\Documents\SAUT';
sujet = {'\001','\002','\003','\004','\005','\006','\007','\008','\009','\010','\011','\012'};
%%
for suj = 1 : 12
    close all
    load([pathname sujet{suj} '\CINE_SAUT.mat'])
    load([pathname sujet{suj} '\EMG_SAUT.mat'])
    
    %%
    
    grpi = 0;
    for grp = 1:4
        grpi = grpi + 1 ;
        mark = find(cellfun(@isempty,strfind({CINE_SAUT(grpi).RawData.Label}, 'meta2_post')') ==0 );
        
        if ~isempty(strfind(CINE_SAUT(grp).Filename, '60'))
            dist = 100;
        end
        if ~isempty(strfind(CINE_SAUT(grp).Filename, '90'))
            dist = 80;
        end 
        
        Cine_meta = -CINE_SAUT(grp).RawData(mark).Data(1,:);
        [pk, lc] = findpeaks(Cine_meta, 'MinPeakDistance', dist);
       
        %% find first peak

        fig = figure(grpi);
        plot(Cine_meta)
        findpeaks(Cine_meta, 'MinPeakDistance', dist)
        dcm_obj = datacursormode(fig);
        set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on')
        waitforbuttonpress
        waitforbuttonpress
        c_info = getCursorInfo(dcm_obj);
        
        [a, b] = min( abs(lc - c_info.Position(1)  ));

        EMG_SAUT(grp).reps = lc(b:2:b+30)/150;
        
        figure(grpi)
        hold on
        line([EMG_SAUT(grp).reps*150 ; EMG_SAUT(grp).reps*150], repmat(ylim,length(EMG_SAUT(grp).reps),1)')
    end  
        save([pathname sujet{suj} '\EMG_SAUT.mat'], 'EMG_SAUT')
    
end
