%% Initialisation

clear all;
close all;
clc;

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))

sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};

pathname='C:\Users\p1218107\Documents\Data_Piano\data_1.5s\'; 

%%


for suj = 1 : 12
    
    load([pathname sujet{suj} '_EMG_GRP.mat'])
    
    for grp = 5 : 8
        for mu = 1 : 10
            
            GRP(grp).EMGred(mu).labels = GRP(grp).EMGal(mu).labels;
            
            GRP(grp).EMGred(mu).data = resample(GRP(grp).EMGal(mu).cycle_ali_resamp',1,4)';

            GRP(grp).EMGred(mu).time = -0.75 + 1.5/788 : 1.5/788 : 0.75;
            
        end
    end
    
    save([pathname sujet{suj} '_EMG_GRP.mat'])
end