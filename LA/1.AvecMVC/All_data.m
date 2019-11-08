%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Normalisation par MVC                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation

clear all;
close all;
clc;

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};

pathname='C:\Users\p1218107\Documents\Data_Piano\LA\'; 

alph = 0.01;


%% Ouverture des fichiers 
load('C:\Users\p1218107\Documents\Data_Piano\MVC\All_MVC01')
load('C:\Users\p1218107\Documents\Data_Piano\MVC\All_MVC1s')
load('C:\Users\p1218107\Documents\Data_Piano\LA\data_newREPS\test_ttt')

%% 
All_MVC1s(1,1) = NaN;
All_MVC1s(7,6) = NaN;

%% Mini 10 cycles
for grp = 1 : 8
    for suj = 1 : 12 
        for mu = 1 : 9
            if size(Test(grp).suj(suj).mu(mu).data_ttt,1)<10
                Fle = Test(grp).suj(suj).mu(mu).data;
                Ext = Test(grp).suj(suj).mu(mu).data_ttt; 
                for i = 1 : size(Ext,1)
                    Fle(Fle(:,1)==Ext(i,1),:) = [];
                end
                
                RMSE = []; 
                for i = 1 : size(Fle,1)
                    RMSE(i) = sqrt( mean( ( Fle(i,:) - mean(Ext) ).^2)); 
                end
            
                [xs, index] = sort(RMSE); 
                Test(grp).suj(suj).mu(mu).data_ttt = [Ext ; Fle(index(1:10-size(Ext,1)) ,: )];
            end
                
        end
    end
end


%% Normalisation + selection
for grp = 1 : 8
    for suj = 1 : 12 
        for mu = 1 : 9
            
            Data = [];
            RMSE = []; 
            
            Data = Test(grp).suj(suj).mu(mu).data_ttt / All_MVC1s(mu,suj);
            
            for i = 1 : size(Data,1)
                RMSE(i) = sqrt( mean( ( Data(i,:) - mean(Data) ).^2)); 
            end
            
            [xs, index] = sort(RMSE); 
            
            Test(grp).suj(suj).mu(mu).data_selec_10_norm = Data(sort(index(1:10))', :);
        end
    end
end

%% Normalisation 
for grp = 1 : 8
    All_data_LA(grp).name = Test(grp).name;
    for mu = 1 : 9
        for suj = 1 : 12   
            All_data_LA(grp).suj(suj).mu(mu).name = Test(grp).suj(suj).mu(mu).muscle;
            All_data_LA(grp).suj(suj).mu(mu).data_norm = Test(grp).suj(suj).mu(mu).data_selec_10_norm; 
            
            All_data_LA(grp).suj_mean(mu).name = Test(grp).suj(suj).mu(mu).muscle;
            All_data_LA(grp).suj_mean(mu).mu(suj,:) = mean(All_data_LA(grp).suj(suj).mu(mu).data_norm);  
        end
    end
end


%%

for suj = 1 : length(sujet)
    for mu = 1 : 9 
        for grp = 1 : 8  
            All_data_LA(grp).suj_mean(mu).mu_res(suj,:) = downsample(All_data_LA(grp).suj_mean(mu).mu(suj,:), 5);
            All_data_LA(grp).suj(suj).mu(mu).data_res =  downsample(All_data_LA(grp).suj(suj).mu(mu).data_norm',5)';
        end
    end 
end

%%

% close all

for mu = 1 : 9
    figure(mu)
    for grp = 1 : 8
        subplot(4,2,grp)
        plot(All_data_LA(grp).suj_mean(mu).mu_res')
        hold on 
    end
    suptitle(All_data_LA(grp).suj_mean(mu).name)
end
%% 

for grp = 1:8
    for suj = 1:12
        for mu = 1 : 9
            All_data_LA(grp).suj(suj).mu(mu).pEMG =  max(All_data_LA(grp).suj(suj).mu(mu).data_res');
            All_data_LA(grp).suj(suj).mu(mu).mEMG =  mean(All_data_LA(grp).suj(suj).mu(mu).data_res');
            All_data_LA(grp).suj_mean(mu).pEMG(suj,1) = nanmean(All_data_LA(grp).suj(suj).mu(mu).pEMG);
            All_data_LA(grp).suj_mean(mu).mEMG(suj,1) = nanmean(All_data_LA(grp).suj(suj).mu(mu).mEMG);
        end
    end
end

%% CI 

% Flexor Extensor 

for grp = 1:8
    for suj = 1:12
        
        Fle = All_data_LA(grp).suj(suj).mu(9).data_res;
        Ext = All_data_LA(grp).suj(suj).mu(8).data_res;
        
        Bic = All_data_LA(grp).suj(suj).mu(2).data_res;
        Tri = All_data_LA(grp).suj(suj).mu(1).data_res;
        
        for i = 1 : 10
            CI_FE(Fle(i,:)<Ext(i,:)) = Fle(i, Fle(i,:)<Ext(i,:));
            CI_FE(Fle(i,:)>Ext(i,:)) = Ext(i, Fle(i,:)>Ext(i,:));            
            All_data_LA(grp).suj(suj).CI_FE(i) = sum(CI_FE)/630;
            
            CI_BT(Bic(i,:)<Tri(i,:)) = Bic(i, Bic(i,:)<Tri(i,:));
            CI_BT(Bic(i,:)>Tri(i,:)) = Tri(i, Bic(i,:)>Tri(i,:)); 
            
            All_data_LA(grp).suj(suj).CI_BT(i) = sum(CI_BT)/630;
        end
        All_data_LA(grp).CI_FE(suj,1) = nanmean(CI_FE);
        All_data_LA(grp).CI_BT(suj,1) = nanmean(CI_BT);
    end
end
