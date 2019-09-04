%% Initialisation

clear all;
close all;
clc;

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))

sujet_name={'001','002','003','004','005','006','007','008','009','010','011','012'};

pathname='C:\Users\p1218107\Documents\Data_Piano\data_1.5s\'; 
cd(pathname)

%% Visualisation all muscles ap ttt sans normalisation

for n = 1 : length(sujet_name)
    
    pathname_emg = [pathname,sujet_name{n},'_EMG_GRP.mat'];
    load(pathname_emg)
    
    for grp = 1 : 8
        for mu = 1 : 10
       
            figure(mu+10*(grp-1))
            %set(gcf, 'Position', get(0, 'Screensize'))
            %suptitle(['grp' num2str(grp) '/mu' num2str(mu)])
            
            subplot(4,3,n)
            plot(GRP(grp).EMGred(mu).data_ttt3')
            legend(num2str(GRP(grp).EMGred(mu).pourc_perte3))
            
        end
    end
    %%
    
    
end

%%
for mu = 1 : 10
    for grp = 1 : 8
        y = [];
        figure(mu+10*(grp-1))
        for n = 1 : length(sujet_name)
            subplot(4,3,n)
             y(n,:)=ylim;
            
        end
        
        ymax = max(y(:,2));
        
        for n = 1 : length(sujet_name) 
            figure(mu+10*(grp-1))
            subplot(4,3,n)
            ylim([0 ymax])
        end
    end
end

%%
for mu = 1 : 10
    for grp = 1 : 8
        figure(mu+10*(grp-1))
        
        suptitle(['grp' num2str(grp) '/mu' num2str(mu)])
        
        set(gcf, 'Position', get(0, 'Screensize'))
        
        saveas(figure(mu+10*(grp-1)), ['C:\Users\p1218107\Documents\Data_Piano\graphiques\muscles\all_muscle all_grps\','grp',num2str(grp),'mu',num2str(mu)],'jpg')

    end 
end
%%

for n = 1 : length(sujet_name)
    
    pathname_emg = [pathname,sujet_name{n},'_EMG_GRP.mat'];
    load(pathname_emg)
    
    for grp = 1 : 8
        for mu = 1 : 10
            
            perte(mu,grp,n) = GRP(grp).EMGred(mu).pourc_perte3;
            
        end
    end
end
perte_bis = perte; 
%%
mu = 8;
grp = 2;
suj = 2;

perte_bis(mu , grp , suj) = 0;

%%
perte_tot2=zeros(10,8);
for n = 1 : 12
    perte_tot2 = perte_tot2 + double(perte_bis(:,:,n)<0.4);
end
            