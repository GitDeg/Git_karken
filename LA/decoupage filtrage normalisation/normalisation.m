%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                Normalisation par moyenne des conditions                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation

clear all;
close all;
clc;

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};

pathname='C:\Users\p1218107\Documents\Data_Piano\data_1.5s\'; 

GrpNames={'Membre supérieur pressé staccato',...
              'Membre supérieur frappé staccato',...
              'Bassin pressé staccato',...
              'Bassin frappé staccato',...
              'Membre supérieur pressé tenue',...
              'Membre supérieur frappé tenue',...
              'Bassin pressé tenue'...
              'Bassin frappé tenue'};

%% Ouverture des fichiers 

for suj = 1 : length(sujet)
    
    pathname_emg = [pathname,sujet{suj},'_EMG_GRP.mat'];
    load(pathname_emg)
    load([pathname 'keep_tot.mat'])
    
    
    %% Normalisation
    for mu = 1 : length(GRP(1).EMGal)
        EMG_mu=[];
        for grp = 1 : length(GRP)
            if keep_tot(mu,grp,suj) == 1
                EMG_mu=[EMG_mu;GRP(grp).EMGred(mu).data_ttt3];
            end
        end
        EMG_mean=nanmean(nanmean(EMG_mu));
%         EMG_sort=sort(EMG_mu(:));
        
        for grp = 1 : length(GRP)
                GRP(grp).EMGred(mu).data_norm3 = [];
                GRP(grp).EMGred(mu).pourc_perte3 = [];
                GRP(grp).EMGred(mu).mean3 = [];
                
            if keep_tot(mu,grp,suj) == 1
                GRP(grp).EMGred(mu).data_norm3 = GRP(grp).EMGred(mu).data_ttt3/EMG_mean;
                GRP(grp).EMGred(mu).pourc_perte3 = size(GRP(grp).EMGred(mu).data_norm3,1)/size(GRP(grp).EMGred(mu).data,1);
                GRP(grp).EMGred(mu).mean3 = mean(GRP(grp).EMGred(mu).data_norm3);
%               GRP(grp).EMGred(mu).max3 = mean(EMG_sort(length(EMG_sort)-round(length(EMG_sort)/1000):end));
            end
        end
    end
    
    save(pathname_emg,'GRP')
end

%% Plot

% ALL_mean_1_5_manu=[];
% ALL_mean_1_5_manu(11).grp1=[];
% ALL_mean_1_5_manu(11).grp2=[];
% ALL_mean_1_5_manu(11).grp3=[];
% ALL_mean_1_5_manu(11).grp4=[];
% ALL_mean_1_5_manu(11).grp5=[];
% ALL_mean_1_5_manu(11).grp6=[];
% ALL_mean_1_5_manu(11).grp7=[];
% ALL_mean_1_5_manu(11).grp8=[];
ALL_grp = [];
ALL_grp.mu = [];
ALL_grp.name = [];

ALL_grp.mu.data = [];
ALL_grp.mu.mean = [];
ALL_grp.mu.name = [];

%%  

GrpNames={'Membre supérieur pressé staccato',...
              'Membre supérieur frappé staccato',...
              'Bassin pressé staccato',...
              'Bassin frappé staccato',...
              'Membre supérieur pressé tenue',...
              'Membre supérieur frappé tenue',...
              'Bassin pressé tenue'...
              'Bassin frappé tenue'};

for suj = 1:12
    pathname_emg = [pathname,sujet{suj},'_EMG_GRP.mat'];
    load(pathname_emg)
%    close all
    %%
    for mu = 1 : 10 
        %%
        for grp = 1 : 8
            
            if suj==1
                
                ALL_grp(grp).mu(mu).data = [GRP(grp).EMGred(mu).data_norm3, suj*ones(size(GRP(grp).EMGred(mu).data_norm3,1),1)];
                ALL_grp(grp).mu(mu).mean = [GRP(grp).EMGred(mu).mean3, suj*ones(size(GRP(grp).EMGred(mu).mean3,1),1)];
                ALL_grp(grp).mu(mu).name = GRP(grp).EMGred(mu).labels;
                ALL_grp(grp).name = GrpNames{grp};
            else
                ALL_grp(grp).mu(mu).data = [ALL_grp(grp).mu(mu).data;...
                                          GRP(grp).EMGred(mu).data_norm3, suj*ones(size(GRP(grp).EMGred(mu).data_norm3,1),1)];
                ALL_grp(grp).mu(mu).mean = [ALL_grp(grp).mu(mu).mean;...
                                          GRP(grp).EMGred(mu).mean3, suj*ones(size(GRP(grp).EMGred(mu).mean3,1),1)];
            end
        end
    end
end
            
            
%             ALL_mean_1_5_manu(1).grp1= GRP(1).Name;
%             ALL_mean_1_5_manu(1).grp2= GRP(2).Name;
%             ALL_mean_1_5_manu(1).grp3= GRP(3).Name;
%             ALL_mean_1_5_manu(1).grp4= GRP(4).Name;
%             ALL_mean_1_5_manu(1).grp5= GRP(5).Name;
%             ALL_mean_1_5_manu(1).grp6= GRP(6).Name;
%             ALL_mean_1_5_manu(1).grp7= GRP(7).Name;
%             ALL_mean_1_5_manu(1).grp8= GRP(8).Name;
%             
%             ALL_mean_1_5_manu(mu+1).grp1= [ALL_mean_1_5_manu(mu+1).grp1;...
%                                   GRP(1).EMGred(mu).mean3, suj*ones(size(GRP(1).EMGred(mu).mean3,1),1)];
%             ALL_mean_1_5_manu(mu+1).grp2= [ALL_mean_1_5_manu(mu+1).grp2;...
%                                   GRP(2).EMGred(mu).mean3, suj*ones(size(GRP(2).EMGred(mu).mean3,1),1)];
%             ALL_mean_1_5_manu(mu+1).grp3= [ALL_mean_1_5_manu(mu+1).grp3;...
%                                   GRP(3).EMGred(mu).mean3, suj*ones(size(GRP(3).EMGred(mu).mean3,1),1)];
%             ALL_mean_1_5_manu(mu+1).grp4= [ALL_mean_1_5_manu(mu+1).grp4;...
%                                   GRP(4).EMGred(mu).mean3, suj*ones(size(GRP(4).EMGred(mu).mean3,1),1)];
%             ALL_mean_1_5_manu(mu+1).grp5= [ALL_mean_1_5_manu(mu+1).grp5;...
%                                   GRP(5).EMGred(mu).mean3, suj*ones(size(GRP(5).EMGred(mu).mean3,1),1)];
%             ALL_mean_1_5_manu(mu+1).grp6= [ALL_mean_1_5_manu(mu+1).grp6;...
%                                   GRP(6).EMGred(mu).mean3, suj*ones(size(GRP(6).EMGred(mu).mean3,1),1)];
%             ALL_mean_1_5_manu(mu+1).grp7= [ALL_mean_1_5_manu(mu+1).grp7;...
%                                   GRP(7).EMGred(mu).mean3, suj*ones(size(GRP(7).EMGred(mu).mean3,1),1)];
%             ALL_mean_1_5_manu(mu+1).grp8= [ALL_mean_1_5_manu(mu+1).grp8;...
%                                   GRP(8).EMGred(mu).mean3, suj*ones(size(GRP(8).EMGred(mu).mean3,1),1)];
%             ALL_mean_1_5_manu(mu+1).Names=GRP(1).EMGal(mu).labels;

       
       
        
%         %%
%         
%     end
%     
% end

%% Plot

for mu = 1 : 10 
    
    figure(1)
    subplot(5,2,mu)
    plot(ALL_mean_1_5_manu(mu + 1).grp1(:,1:end-1)')
    
    figure(2)
    subplot(5,2,mu)
    plot(ALL_mean_1_5_manu(mu + 1).grp2(:,1:end-1)')
    
    figure(3)
    subplot(5,2,mu)
    plot(ALL_mean_1_5_manu(mu + 1).grp3(:,1:end-1)')
    
    figure(4)
    subplot(5,2,mu)
    plot(ALL_mean_1_5_manu(mu + 1).grp4(:,1:end-1)')
    
    figure(5)
    subplot(5,2,mu)
    plot(ALL_mean_1_5_manu(mu + 1).grp5(:,1:end-1)')
    
    figure(6)
    subplot(5,2,mu)
    plot(ALL_mean_1_5_manu(mu + 1).grp6(:,1:end-1)')
    
    figure(7)
    subplot(5,2,mu)
    plot(ALL_mean_1_5_manu(mu + 1).grp7(:,1:end-1)')
    
    figure(8)
    subplot(5,2,mu)
    plot(ALL_mean_1_5_manu(mu + 1).grp8(:,1:end-1)')
end
    



% %% retrait data 
% % Sujet 1 grp4 grand dent(10)
% % Sujet 2 grp2 grand pec(8)
% % Sujet 9 grp4 extenseur (3) et deltmed(7)
%% 
grp = 2;
mu = 8;

for suj = 2
%     
%     pathname_emg = [pathname,sujet{n},'\EMG_grp_30ms.mat'];
%     load(pathname_emg)
    
    GRP(grp).EMGred(mu).data_norm=NaN(1,788);
%     GRP(grp).EMGal(mu).mean_norm=NaN(1,1049);
%     GRP(grp).EMGal(mu).cycl_norm=NaN(1,1049);
%     
%    save(pathname_emg,'GRP')
end
%     
%     
% 
%ALL_mean(mu+1).grp4(any(isnan(ALL_mean(mu+1).grp4),2),:) = nanmean(ALL_mean(mu+1).grp4);