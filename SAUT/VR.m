%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                           VR_intra SAUT                           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc 

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
%%
% pathname = '\\10.89.24.15\j\Valentin\SAUT';
pathname = 'C:\Users\p1218107\Documents\SAUT';
sujet = {'\001','\002','\003','\004','\005','\006','\007','\008','\009','\010','\011','\012'};

load([pathname '\All_data.mat'])

%% Sélection des Datas à traiter  
%%%%%%%%%%%%% EMG %%%%%%%%%%%%%%

for suj = 1 : 12 
    for grp = 1 : 6
        for mu = 1 : 9    
            %% Variance ratio (coefficient de dispertion) INTRA cycle : entre les cycles d'un meme muscle 
            % On utilise toutes les données pour calculer la dispertion des cycles
            % autour de la moyenne. ce coefficient permet de comparer des variations
            % lorsque les moyennes sont differentes. 

            
            Xij = All_data(grp).EMG(suj).mu(mu).data  ;
            nb_frame = length(Xij); % nb de frame
            n_data(grp,mu,suj) = size(Xij,1);

            Xi = nanmean(Xij);
            X = nanmean(nanmean(Xi));

            num = nansum( nansum( (Xij - Xi).^2  / (nb_frame*(n_data(grp,mu,suj)-1))));
            den = nansum( nansum( (Xij - X).^2 / (nb_frame*n_data(grp,mu,suj) -1)));

            All_data(grp).EMG(suj).mu(mu).VR_intra = num / den ;% nb de cycle 
            
            
        end
    end
end

%% 

for grp = 1 : 6
    for mu = 1 : 9   
        %% Variance ratio (coefficient de dispertion) INTER individuel : entre les sujets  
        % Cette fois on utilise les données moyennes pour calculer la variation
        % entre les sujet réalisant un meme test. 

        Xij = All_data(grp).EMG(1).mean_EMG(mu).data      ;
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
            % On utilise toutes les données pour calculer la dispertion des cycles
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
            
            %%
            Yij = All_data(grp).CINE(suj).mark(mark).data((1:10)*3-1,:)    ;
            nb_frame = length(Yij); % nb de frame
            n_data = size(Yij,1);

            Yi = nanmean(Yij);
            Y = nanmean(nanmean(Yi));

            num_Y = nansum( nansum( (Yij - Yi).^2  / (nb_frame*(n_data-1))));
            den_Y = nansum( nansum( (Yij - Y).^2 / (nb_frame*n_data -1)));

            All_data(grp).CINE(suj).mark(mark).VR_intra_Y = num_Y / den_Y ;% nb de cycle 
            
            %%
            Zij = All_data(grp).CINE(suj).mark(mark).data((1:10)*3,:)    ;
            nb_frame = length(Zij); % nb de frame
            n_data = size(Zij,1);

            Zi = nanmean(Zij);
            Z = nanmean(nanmean(Zi));

            num_Z = nansum( nansum( (Zij - Zi).^2  / (nb_frame*(n_data -1))));
            den_Z = nansum( nansum( (Zij - Z).^2 / (nb_frame*n_data -1)));

            All_data(grp).CINE(suj).mark(mark).VR_intra_Z = num_Z / den_Z ;% nb de cycle 
            
            
        end
    end
end

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

subplot(3,2,2)
boxplot(boxplot_data_EMG(ordre_mu,:,2)',{All_data(1).EMG(1).mu(ordre_mu).name})
title(All_data(2).Filename)

subplot(3,2,3)
boxplot(boxplot_data_EMG(ordre_mu,:,3)',{All_data(1).EMG(1).mu(ordre_mu).name})
title(All_data(3).Filename)

subplot(3,2,4)
boxplot(boxplot_data_EMG(ordre_mu,:,4)',{All_data(1).EMG(1).mu(ordre_mu).name})
title(All_data(4).Filename)

subplot(3,2,5)
boxplot(boxplot_data_EMG(ordre_mu,:,5)',{All_data(1).EMG(1).mu(ordre_mu).name})
title(All_data(5).Filename)

subplot(3,2,6)
boxplot(boxplot_data_EMG(ordre_mu,:,6)',{All_data(1).EMG(1).mu(ordre_mu).name})
title(All_data(6).Filename)

%%
figure 
mui = 0;
for mu = ordre_mu 
    mui = mui + 1;
    subplot(3,3,mui)
    boxplot(reshape(boxplot_data_EMG(mu,:,:),12,6))
    title({All_data(1).EMG(1).mu(mu).name})
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
title('X antero-posterieur')
subplot(3,1,2)
boxplot([reshape(boxplot_data_META5(2,:,:),12,4), reshape(boxplot_data_META2(2,1:12,5:6),12,2)])
title('Y medio-lateral')
subplot(3,1,3)
boxplot([reshape(boxplot_data_META5(3,:,:),12,4), reshape(boxplot_data_META2(3,1:12,5:6),12,2)])
title('Z cranio-plantaire')

%%
close all
figure

grp = 5;
for suj = 1 : 12
    figure(suj) ; hold on; 

%    plot3(All_data(1).CINE(suj).mark(1).mean_data(1,:), All_data(1).CINE(suj).mark(1).mean_data(2,:), All_data(1).CINE(suj).mark(1).mean_data(3,:), 'color', 'r', 'linewidth' , 2);
    for cycle = 1 : 10
        plot3(All_data(grp).CINE(suj).mark(1).data(cycle*3-2,:), All_data(grp).CINE(suj).mark(1).data(cycle*3-1,:), All_data(grp).CINE(suj).mark(1).data(cycle*3,:), 'color', 'b');
    end
    grid on
    axis equal
end

%%

