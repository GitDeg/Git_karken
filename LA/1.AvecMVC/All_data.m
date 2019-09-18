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
load('C:\Users\p1218107\Documents\Data_Piano\MVC\All_MVC')
load('C:\Users\p1218107\Documents\Data_Piano\LA\data_newREPS\test_ttt')

%% Mini 10 cycles
for grp = 1 : 8
    for suj = 1 : 12 
        for mu = 1 : 9
            if size(Test(grp).suj(suj).mu(mu).data_ttt,1)<10
                A = Test(grp).suj(suj).mu(mu).data;
                B = Test(grp).suj(suj).mu(mu).data_ttt; 
                for i = 1 : size(B,1)
                    A(A(:,1)==B(i,1),:) = [];
                end
                
                RMSE = []; 
                for i = 1 : size(A,1)
                    RMSE(i) = sqrt( mean( ( A(i,:) - mean(B) ).^2)); 
                end
            
                [xs, index] = sort(RMSE); 
                Test(grp).suj(suj).mu(mu).data_ttt = [B ; A(index(1:10-size(B,1)) ,: )];
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
            
            Data = Test(grp).suj(suj).mu(mu).data_ttt / All_MVC(mu,suj);
            
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

close all

for mu = 1 : 9
    figure(mu)
    for grp = 1 : 8
        subplot(4,2,grp)
        plot(All_data_LA(grp).suj_mean(mu).mu_res')
        hold on 
    end
    suptitle(All_data_LA(grp).suj_mean(mu).name)
end
