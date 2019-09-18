%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                           VR_intra SAUT                           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc 

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
%%
% pathname = '\\10.89.24.15\j\Valentin\SAUT';
pathname = 'C:\Users\p1218107\Documents\Data_Piano\MVC';
sujet = {'\001','\002','\003','\004','\005','\006','\007','\008','\009','\010','\011','\012'};

%%

suj = 2;% : 12
    
    load([pathname sujet{suj} '\MVC.mat'])
    EMG = [];
    for grp = 1 : length(MVC)
        EMG = [EMG , MVC(grp).RawData];
    end
    
    figure(suj)
    for i = 1 : 10
        subplot(5,2,i)
        plot(EMG(i,:))
        title(MVC(1).Labels(i))
    end
    
    
    %%
    mu_selec = [2 3];
    
    for mu = 1:10
        MVC_selec(mu).name = MVC(1).Labels{mu};
        
        if sum(mu == mu_selec)~=0
            
            temp = sort( abs(EMG(mu,:)));
            
            MVC_selec(mu).MVC = mean(temp( round((length(temp)-0.001*length(temp))) : length(temp)));
        end
        
    end
    
    %%
    mu_incom = [1 4 5 6 7 8 9 10];
    for mu = 10
        if sum(mu == mu_incom)~=0
            interv = [1:14800 15500:length(EMG(mu,:))];
            temp = sort( abs(EMG(mu, interv)));
            MVC_selec(mu).MVC = mean(temp( round((length(temp)-0.001*length(temp))) : length(temp)));
        end
    end
    
    
    %%
    
save([pathname sujet{suj} '\MVC_selec.mat'], 'MVC_selec')

%%

for suj = 1 : 12
    load([pathname sujet{suj} '\MVC_selec.mat'])
    
    for i = 1 : 10
        if isempty(MVC_selec(i).MVC)
            
            MVC_selec(i).MVC = NaN;
            
        end
    end
    
    All_MVC(1:8,suj) = [MVC_selec(1:8).MVC]';
    All_MVC(9,suj) = nanmean([MVC_selec(9:10).MVC]);
end

save([pathname '\All_MVC.mat'], 'All_MVC')

%%

nanmean(All_MVC,2)