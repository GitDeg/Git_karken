%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                           VR_intra SAUT                           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc 

%% Ubuntu 

% addpath(genpath('/media/valentin/060C6AFC0C6AE5E1/Users/p1218107/Documents/MATLAB'))
% 
% pathname = '/media/valentin/060C6AFC0C6AE5E1/Users/p1218107/Documents/Data_Piano/SAUT';
% sujet = {'/001','/002','/003','/004','/005','/006','/007','/008','/009','/010','/011','/012'};
% 
% load([pathname '/All_data1s.mat'])

%% Windows
% pathname = '\\10.89.24.15\j\Valentin\SAUT';
addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
pathname = 'C:\Users\p1218107\Documents\Data_Piano\SAUT';

sujet = {'\001','\002','\003','\004','\005','\006','\007','\008','\009','\010','\011','\012'};

load([pathname '\All_data1s_xy.mat'])

%% S�lection des Datas � traiter  
%%%%%%%%%%%%% EMG %%%%%%%%%%%%%%

for suj = 1 : 12 
    for grp = 1 : 6
        for mu = 1 : 9    
            %% Variance ratio (coefficient de dispertion) INTRA cycle : entre les cycles d'un meme muscle 
            % On utilise toutes les donn�es pour calculer la dispertion des cycles
            % autour de la moyenne. ce coefficient permet de comparer des variations
            % lorsque les moyennes sont differentes. 

            
            Xij = All_data(grp).EMG(suj).mu(mu).data_MVC  ;
            nb_frame = length(Xij); % nb de frame
            n_data(grp,mu,suj) = size(Xij,1);

            Xi = nanmean(Xij);
            X = nanmean(nanmean(Xi));

            num = nansum( nansum( (Xij - Xi).^2  / (nb_frame*(n_data(grp,mu,suj)-1))));
            den = nansum( nansum( (Xij - X).^2 / (nb_frame*n_data(grp,mu,suj) -1)));

%             All_data(grp).EMG(suj).mu(mu).VR_intra = num / den ;% nb de cycle 
%             
%             All_data(grp).EMG(suj).mu(mu).Peak = max(Xij');
%             
            %% CV
            All_data(grp).EMG(suj).mu(mu).CV = sqrt(nanmean( nanstd(Xij).*nanstd(Xij) ))/nanmean(nanmean(Xij));
            All_data(grp).EMG(1).mean_EMG(mu).CVi(suj,:) = nanstd(Xij)./nanmean(Xij);
%             
%             %% CQV
%             All_data(grp).EMG(suj).mu(mu).CQV = (prctile(Xij,75) - prctile(Xij,25)) ./ (prctile(Xij,75) + prctile(Xij,25));
%             
        end
    end
end

%%
grp = 5
for mu = 1 : 9
    for suj = 1 : 12
        peak_grp5(mu).mu(:,suj) =  All_data(grp).EMG(suj).mu(mu).Peak;

    end
end


for mu = 1 : 9
    figure(mu)
    boxplot([peak_grp1(mu).mu(:),peak_grp5(mu).mu(:)])
end
% figure
% subplot(3,1,1); plot(All_data(1).EMG(1).mu(1).CVi);% hold on; plot(mean(All_data(1).EMG(1).mu(1).data)); hold on; plot(std(All_data(1).EMG(1).mu(1).data))
% subplot(3,1,2); plot(All_data(1).EMG(1).mu(4).CVi);% hold on; plot(mean(All_data(1).EMG(1).mu(8).data)); hold on; plot(std(All_data(1).EMG(1).mu(8).data))
% subplot(3,1,3); plot(All_data(1).EMG(1).mu(9).CVi);% hold on; plot(mean(All_data(1).EMG(1).mu(9).data)); hold on; plot(std(All_data(1).EMG(1).mu(9).data))
% figure 
% subplot(3,1,1); plot(All_data(5).EMG(1).mu(1).CVi);% hold on; plot(mean(All_data(5).EMG(1).mu(1).data)); hold on; plot(std(All_data(5).EMG(1).mu(1).data))
% subplot(3,1,2); plot(All_data(5).EMG(1).mu(4).CVi);% hold on; plot(mean(All_data(5).EMG(1).mu(8).data)); hold on; plot(std(All_data(5).EMG(1).mu(8).data))
% subplot(3,1,3); plot(All_data(5).EMG(1).mu(9).CVi);% hold on; plot(mean(All_data(5).EMG(1).mu(9).data)); hold on; plot(std(All_data(5).EMG(1).mu(9).data))
% 
% 
% figure
% subplot(3,1,1); plot(All_data(1).EMG(1).mu(1).CQV);% hold on; plot(mean(All_data(1).EMG(1).mu(1).data)); hold on; plot(std(All_data(1).EMG(1).mu(1).data))
% subplot(3,1,2); plot(All_data(1).EMG(1).mu(4).CQV);% hold on; plot(mean(All_data(1).EMG(1).mu(8).data)); hold on; plot(std(All_data(1).EMG(1).mu(8).data))
% subplot(3,1,3); plot(All_data(1).EMG(1).mu(9).CQV);% hold on; plot(mean(All_data(1).EMG(1).mu(9).data)); hold on; plot(std(All_data(1).EMG(1).mu(9).data))
% figure 
% subplot(3,1,1); plot(All_data(5).EMG(1).mu(1).CQV);% hold on; plot(mean(All_data(5).EMG(1).mu(1).data)); hold on; plot(std(All_data(5).EMG(1).mu(1).data))
% subplot(3,1,2); plot(All_data(5).EMG(1).mu(4).CQV);% hold on; plot(mean(All_data(5).EMG(1).mu(8).data)); hold on; plot(std(All_data(5).EMG(1).mu(8).data))
% subplot(3,1,3); plot(All_data(5).EMG(1).mu(9).CQV);% hold on; plot(mean(All_data(5).EMG(1).mu(9).data)); hold on; plot(std(All_data(5).EMG(1).mu(9).data))
% 

%% 

for grp = 1 : 6
    for mu = 1 : 9   
        %% Variance ratio (coefficient de dispertion) INTER individuel : entre les sujets  
        % Cette fois on utilise les donn�es moyennes pour calculer la variation
        % entre les sujet r�alisant un meme test. 

        Xij = All_data(grp).EMG(1).mean_EMG(mu).data_MVC      ;
        Xi = nanmean(Xij);
        X = nanmean(Xi);

        %n_mean = size(ALL_grp(grp).mu(mu).mean(:,1:id_suj),1);
        n_mean = size(Xij,1); % nb de sujet 
        nb_frame = size(Xij,2); % nb de frames

        num = nansum( nansum( (Xij - Xi).^2)) / (nb_frame*(n_mean-1));
        den = nansum( nansum( (Xij - X).^2))/ (nb_frame*n_mean -1);

        All_data(grp).EMG(1).mean_EMG(mu).VR_inter = num / den ;
    end
end

%% 
%%%%%%%%% CINE %%%%%%%%%

for grp = 1 : length(All_data) 
    for suj = 1 : length(All_data(grp).CINE)
        for mark = 1 : length(All_data(grp).CINE(suj).mark  )   
            %% Variance ratio (coefficient de dispertion) INTRA cycle : entre les cycles d'un meme muscle 
            % On utilise toutes les donn�es pour calculer la dispertion des cycles
            % autour de la moyenne. ce coefficient permet de comparer des variations
            % lorsque les moyennes sont differentes. 

            
            Xij = All_data(grp).CINE(suj).mark(mark).data((1:10)*3-2,:)    ;
            nb_frame = length(Xij); % nb de frame
            n_data = size(Xij,1);

            Xi = nanmean(Xij);
            X = nanmean(nanmean(Xi));

            num_X = nansum( nansum( (Xij - Xi).^2  / (nb_frame*(n_data-1))));
            den_X = nansum( nansum( (Xij - X).^2 / (nb_frame*n_data -1)));

            All_data(grp).CINE(suj).mark(mark).VR_intra_X = num_X / den_X ;% nb de cycle 
%             
            All_data(grp).CINE(suj).mark(mark).CV_X = sqrt(nanmean( nanstd(Xij).*nanstd(Xij) ))/nanmean(nanmean(abs(Xij)));
            All_data(grp).CINE(suj).mark(mark).CVi_X = nanstd(Xij)./nanmean(abs(Xij));
            All_data(grp).CINE(suj).mark(mark).CQV_X = (prctile(Xij,75) - prctile(Xij,25)) ./ (prctile(Xij,75) + prctile(Xij,25));
            
            %%
            Yij = All_data(grp).CINE(suj).mark(mark).data((1:10)*3-1,:)    ;
            nb_frame = length(Yij); % nb de frame
            n_data = size(Yij,1);

            Yi = nanmean(Yij);
            Y = nanmean(nanmean(Yi));

            num_Y = nansum( nansum( (Yij - Yi).^2  / (nb_frame*(n_data-1))));
            den_Y = nansum( nansum( (Yij - Y).^2 / (nb_frame*n_data -1)));

            All_data(grp).CINE(suj).mark(mark).VR_intra_Y = num_Y / den_Y ;% nb de cycle 
%             
            All_data(grp).CINE(suj).mark(mark).CV_Y = sqrt(nanmean( nanstd(Yij).*nanstd(Yij) ))/nanmean(nanmean(abs(Yij)));
            All_data(grp).CINE(suj).mark(mark).CVi_Y = nanstd(Yij)./nanmean(abs(Yij));
            All_data(grp).CINE(suj).mark(mark).CQV_Y = (prctile(Yij,75) - prctile(Yij,25)) ./ (prctile(Yij,75) + prctile(Yij,25));
            
            %%
            Zij = All_data(grp).CINE(suj).mark(mark).data((1:10)*3,:)    ;
            nb_frame = length(Zij); % nb de frame
            n_data = size(Zij,1);

            Zi = nanmean(Zij);
            Z = nanmean(nanmean(Zi));

            num_Z = nansum( nansum( (Zij - Zi).^2  / (nb_frame*(n_data -1))));
            den_Z = nansum( nansum( (Zij - Z).^2 / (nb_frame*n_data -1)));

            All_data(grp).CINE(suj).mark(mark).VR_intra_Z = num_Z / den_Z ;% nb de cycle 
%             
            All_data(grp).CINE(suj).mark(mark).CV_Z = sqrt(nanmean( nanstd(Zij).*nanstd(Zij) ))/nanmean(nanmean(abs(Zij)));
            All_data(grp).CINE(suj).mark(mark).CVi_Z = nanstd(Zij)./nanmean(abs(Zij));
            All_data(grp).CINE(suj).mark(mark).CQV_Z = (prctile(Zij,75) - prctile(Zij,25)) ./ (prctile(Zij,75) + prctile(Zij,25));
            
%             
        end
    end
end

%%
close all
figure(1)
subplot(3,1,1); plot(All_data(1).CINE(1).mark(1).CVi_X);
subplot(3,1,2); plot(All_data(1).CINE(1).mark(1).CVi_Y); 
subplot(3,1,3); plot(All_data(1).CINE(1).mark(1).CVi_Z);
figure(2) 
subplot(3,1,1); plot(All_data(3).CINE(1).mark(1).CVi_X); 
subplot(3,1,2); plot(All_data(3).CINE(1).mark(1).CVi_Y); 
subplot(3,1,3); plot(All_data(3).CINE(1).mark(1).CVi_Z);
figure(3) 
subplot(3,1,1); plot(All_data(5).CINE(1).mark(1).CVi_X); 
subplot(3,1,2); plot(All_data(5).CINE(1).mark(1).CVi_Y); 
subplot(3,1,3); plot(All_data(5).CINE(1).mark(1).CVi_Z);
%%
close all
figure
subplot(3,1,1); plot(All_data(1).CINE(1).mark(1).CQV_X);
subplot(3,1,2); plot(All_data(1).CINE(1).mark(1).CQV_Y); 
subplot(3,1,3); plot(All_data(1).CINE(1).mark(1).CQV_Z); 
figure
subplot(3,1,1); plot(All_data(3).CINE(1).mark(1).CQV_X);
subplot(3,1,2); plot(All_data(3).CINE(1).mark(1).CQV_Y); 
subplot(3,1,3); plot(All_data(3).CINE(1).mark(1).CQV_Z); 
figure 
subplot(3,1,1); plot(All_data(5).CINE(1).mark(1).CQV_X); 
subplot(3,1,2); plot(All_data(5).CINE(1).mark(1).CQV_Y); 
subplot(3,1,3); plot(All_data(5).CINE(1).mark(1).CQV_Z);
%%
for suj = 1 : 12
    for grp = 1 : 6
            
            boxplot_data_EMG(:,suj,grp) = [All_data(grp).EMG(suj).mu.VR_intra]';
            
    end
end

%%
ordre_mu = [9 8 2 1 3 4 7 5 6];
figure 
subplot(3,2,1)
boxplot(boxplot_data_EMG(ordre_mu,:,1)',{All_data(1).EMG(1).mu(ordre_mu).name})
title(All_data(1).Filename)
ylim([0 1])

subplot(3,2,2)
boxplot(boxplot_data_EMG(ordre_mu,:,2)',{All_data(1).EMG(1).mu(ordre_mu).name})
title(All_data(2).Filename)
ylim([0 1])

subplot(3,2,3)
boxplot(boxplot_data_EMG(ordre_mu,:,3)',{All_data(1).EMG(1).mu(ordre_mu).name})
title(All_data(3).Filename)
ylim([0 1])

subplot(3,2,4)
boxplot(boxplot_data_EMG(ordre_mu,:,4)',{All_data(1).EMG(1).mu(ordre_mu).name})
title(All_data(4).Filename)
ylim([0 1])

subplot(3,2,5)
boxplot(boxplot_data_EMG(ordre_mu,:,5)',{All_data(1).EMG(1).mu(ordre_mu).name})
title(All_data(5).Filename)
ylim([0 1])

subplot(3,2,6)
boxplot(boxplot_data_EMG(ordre_mu,:,6)',{All_data(1).EMG(1).mu(ordre_mu).name})
title(All_data(6).Filename)
ylim([0 1])

%%
figure 
mui = 0;
for mu = ordre_mu 
    mui = mui + 1;
    subplot(3,3,mui)
    boxplot(reshape(boxplot_data_EMG(mu,:,:),12,6))
    title({All_data(1).EMG(1).mu(mu).name})
    ylim([0 1])
end


%%
for suj = 1 : 12
    for grp = 1 : 4
       
        boxplot_data_META1(1,suj,grp) = All_data(grp).CINE(suj).mark(1).VR_intra_X;
        boxplot_data_META1(2,suj,grp) = All_data(grp).CINE(suj).mark(1).VR_intra_Y;
        boxplot_data_META1(3,suj,grp) = All_data(grp).CINE(suj).mark(1).VR_intra_Z;
        
        boxplot_data_META5(1,suj,grp) = All_data(grp).CINE(suj).mark(2).VR_intra_X;
        boxplot_data_META5(2,suj,grp) = All_data(grp).CINE(suj).mark(2).VR_intra_Y;
        boxplot_data_META5(3,suj,grp) = All_data(grp).CINE(suj).mark(2).VR_intra_Z;
        
    end
end

for suj = 1 : 12
    for grp = 5:6
        boxplot_data_META2(1,suj,grp) = All_data(grp).CINE(suj).mark(1).VR_intra_X;
        boxplot_data_META2(2,suj,grp) = All_data(grp).CINE(suj).mark(1).VR_intra_Y;
        boxplot_data_META2(3,suj,grp) = All_data(grp).CINE(suj).mark(1).VR_intra_Z;
    end
end

%%

figure 
subplot(3,2,1)
boxplot(boxplot_data_META5(:,:,1)',{'X' 'Y' 'Z'})
ylim([0 1])
title(All_data(1).Filename)

subplot(3,2,2)
boxplot(boxplot_data_META5(:,:,2)',{'X' 'Y' 'Z'})
ylim([0 1])
title(All_data(2).Filename)

subplot(3,2,3)
boxplot(boxplot_data_META5(:,:,3)',{'X' 'Y' 'Z'})
ylim([0 1])
title(All_data(3).Filename)

subplot(3,2,4)
boxplot(boxplot_data_META5(:,:,4)',{'X' 'Y' 'Z'})
ylim([0 1])
title(All_data(4).Filename)

subplot(3,2,5)
boxplot(boxplot_data_META2(:,:,5)',{'X' 'Y' 'Z'})
ylim([0 1])
title(All_data(5).Filename)

subplot(3,2,6)
boxplot(boxplot_data_META2(:,:,6)',{'X' 'Y' 'Z'})
ylim([0 1])
title(All_data(6).Filename)

%%
figure
hold on
subplot(3,1,1)
boxplot([reshape(boxplot_data_META5(1,:,:),12,4), reshape(boxplot_data_META2(1,1:12,5:6),12,2)])
title('X medio-lateral')
ylim([0 1])
subplot(3,1,2)
boxplot([reshape(boxplot_data_META5(2,:,:),12,4), reshape(boxplot_data_META2(2,1:12,5:6),12,2)])
title('Y antero-posterieur')
ylim([0 1])
subplot(3,1,3)
boxplot([reshape(boxplot_data_META5(3,:,:),12,4), reshape(boxplot_data_META2(3,1:12,5:6),12,2)])
title('Z cranio-plantaire')
ylim([0 1])

%%

close all
c = {'r', 'b', 'g'};
ind = 1;
for grp =  [1 4 6]
    for suj = 11% : 12
        figure(suj)
        %subplot(3,2,grp);
        if grp < 5
            mark = 5;
        else
            mark = 1;
        end
        hold on; 
        plot3(All_data(grp).CINE(suj).mark(mark).mean_data(1,:), All_data(grp).CINE(suj).mark(mark).mean_data(2,:), All_data(grp).CINE(suj).mark(mark).mean_data(3,:), 'color', c{ind}, 'linewidth' , 2);
    %    for cycle = 1 : 10
    %        plot3(All_data(grp).CINE(suj).mark(1).data(cycle*3-2,:), All_data(grp).CINE(suj).mark(1).data(cycle*3-1,:), All_data(grp).CINE(suj).mark(1).data(cycle*3,:), 'color', 'b');
    %    end
        grid on
        axis equal
        ind = ind + 1;
    end
end

%%

