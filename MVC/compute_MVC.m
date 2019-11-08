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

for suj = 1 : 12
    close all
    
    load([pathname sujet{suj} '\MVC.mat'])
    EMG = [];
    for grp = 1 : length(MVC)
        EMG = [EMG , MVC(grp).RawData];
    end

    for mu = 1 : 10
        figure(mu)
        MVC_corr(suj).mu(mu).data = EMG(mu,:);
        plot(MVC_corr(suj).mu(mu).data)
        title([MVC(1).Labels(mu), 'sujet ', sujet{suj}])

        sat = 1;
        sat = input('saturation ? 1 = oui / 0 = non ');
        
        while sat == 1 

            s = imrect;
            P = getPosition(s);
            if (round(P(1)) +round(P(3))) > length(MVC_corr(suj).mu(mu).data)
                P(3) = length(MVC_corr(suj).mu(mu).data) - round(P(1));
            end
            MVC_corr(suj).mu(mu).data(round(P(1)): round(P(1)) + round(P(3))) = [];
            figure(mu)
            plot(MVC_corr(suj).mu(mu).data)
            title([MVC(1).Labels(mu), 'sujet ', sujet{suj}])
            
            sat = input('saturation ? 1 = oui / 0 = non ');
        end
    end
end


    
    %%
    
load([pathname '\MVC_corr.mat'])

F = 2250;

for suj = 1 : 12
    
    for mu = 1:8
        
        A = sort(MVC_corr(suj).mu(mu).data);
        n = length(A);
        
        All_MVC1s(mu,suj) = nanmedian(A(n-2250 : n));
    end
    
    A1 = sort(MVC_corr(suj).mu(9).data);
    A2 = sort(MVC_corr(suj).mu(10).data);
    n1 = length(A1);
    n2 = length(A2);
    
    All_MVC1s(9,suj) = nanmean( [nanmedian(A1(n1-2250 : n1)), nanmedian(A2(n2-2250 : n2)) ] );
end


%% 
All_MVC1s(3,1) = NaN;
All_MVC1s(7,2) = NaN;
save([pathname '\All_MVC1s.mat'], 'All_MVC1s')

%%