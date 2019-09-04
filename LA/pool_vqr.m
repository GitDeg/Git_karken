%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                Normalisation par moyenne des conditions                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation

clear all;
close all;
clc;

%addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB\spm1dmatlab-master'))
sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};

pathname='C:\Users\p1218107\Documents\Data_Piano\'; 




%% Plot

ALL_pool_var=[];
ALL_pool_var(13).grp1=[];
ALL_pool_var(13).grp2=[];
ALL_pool_var(13).grp3=[];
ALL_pool_var(13).grp4=[];


   

%%  

for n = 1 : 12
    pathname_emg = [pathname,sujet{n},'\EMG_grp_30ms.mat'];
    load(pathname_emg)
    close all
    %%

    A = [];
    A.grp1 = [];
    A.grp2 = []; 
    A.grp3 = [];
    A.grp4 = []; 
    
    Sz = [];
    Sz.grp1 = 0;
    Sz.grp2 = 0;
    Sz.grp3 = 0;
    Sz.grp4 = 0;
    for mu = 1 : length(GRP(1).EMGal)  
        
            ALL_pool_var(1).grp1= GRP(1).Name;
            ALL_pool_var(1).grp2= GRP(2).Name;
            ALL_pool_var(1).grp3= GRP(3).Name;
            ALL_pool_var(1).grp4= GRP(4).Name;

            
            A.grp1= [A.grp1;...
                                  (size(GRP(1).EMGal(mu).cycle_ali_ttt,1)-1)*var(GRP(1).EMGal(mu).cycle_ali_ttt,0,1)./...
                                  mean(GRP(1).EMGal(mu).cycle_ali_ttt,1)];
            A.grp2= [A.grp2;...
                                  (size(GRP(2).EMGal(mu).cycle_ali_ttt,1)-1)*var(GRP(2).EMGal(mu).cycle_ali_ttt,0,1)./...
                                  mean(GRP(2).EMGal(mu).cycle_ali_ttt,1)];
            A.grp3= [A.grp3;...
                                  (size(GRP(3).EMGal(mu).cycle_ali_ttt,1)-1)*var(GRP(3).EMGal(mu).cycle_ali_ttt,0,1)./...
                                  mean(GRP(3).EMGal(mu).cycle_ali_ttt,1)];
            A.grp4= [A.grp4;...
                                  (size(GRP(4).EMGal(mu).cycle_ali_ttt,1)-1)*var(GRP(4).EMGal(mu).cycle_ali_ttt,0,1)./...
                                  mean(GRP(4).EMGal(mu).cycle_ali_ttt,1)];
            
            Sz.grp1 = Sz.grp1 + size(GRP(1).EMGal(mu).cycle_ali_ttt,1);
            Sz.grp2 = Sz.grp2 + size(GRP(2).EMGal(mu).cycle_ali_ttt,1);
            Sz.grp3 = Sz.grp3 + size(GRP(3).EMGal(mu).cycle_ali_ttt,1);
            Sz.grp4 = Sz.grp4 + size(GRP(4).EMGal(mu).cycle_ali_ttt,1);

    end
    
            ALL_pool_var(n+1).grp1= nansum(A.grp1)/(Sz.grp1 - 10);
            ALL_pool_var(n+1).grp2= nansum(A.grp2)/(Sz.grp2 - 10);
            ALL_pool_var(n+1).grp3= nansum(A.grp3)/(Sz.grp3 - 10);
            ALL_pool_var(n+1).grp4= nansum(A.grp4)/(Sz.grp4 - 10);

end
 

save('C:\Users\p1218107\Documents\Data_Piano\All_pool_var.mat','ALL_pool_var')
% ALL_mean(mu+1).grp4(any(isnan(ALL_mean(mu+1).grp4),2),:) = nanmean(ALL_mean(mu+1).grp4);