%% Clustering par minimisation CV_inter ou par minimisation de rmse
clear all;
close all;
clc;

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))

load('C:\Users\p1218107\Documents\Data_Piano\data_newREPS\All_data2.mat')

%% Variance ratio (coefficient de dispertion) INTER individuel : entre les sujets  
        % Cette fois on utilise les données moyennes pour calculer la variation
        % entre les sujet réalisant un meme test. 

for grp = 1 %: 8
    for mu = 1 : 9
        Xij = All_data(grp).suj_mean(mu).mu  ;
        Xi = nanmean(Xij);
        X = nanmean(Xi);

        n_mean = size(Xij,1); % nb de sujet 
        nb_frame = size(Xij,2); % nb de frames

        num = nansum( nansum( (Xij - Xi).^2  / (nb_frame*(n_mean-1))));
        den = nansum( nansum( (Xij - X).^2 / (nb_frame*n_mean -1)));

        All_data(grp).suj_mean(mu).VR_inter = num / den ;
    end
end

%% rmse + corr

for grp = 1 : 8
    for mu = 1 : 9
            DATA = All_data(grp).suj_mean(mu).mu'; 
            
        for i = 1 : size(DATA,2)
            for j = 1 : size(DATA,2)
                All_data(grp).suj_mean(mu).rmse(i,j) = sqrt( sum( (DATA(:,i) - DATA(:,j)).^2))/numel(DATA(:,i));
            end
        end       
         All_data(grp).suj_mean(mu).corr = corrcoef(DATA);
    end
end



%% Minimisation CV_inter

for grp = 1 : 8
    for mu = 1 : 9
        
        suj_set = 1:12;
        VR_inter = [];
        C = nchoosek(suj_set, 6);
        
        for i = 1:length(C)
            Xij = All_data(grp).suj_mean(mu).mu(C(i,:),:)  ;
            Xi = nanmean(Xij);
            X = nanmean(Xi);

            n_mean = size(Xij,1); % nb de sujet 
            nb_frame = size(Xij,2); % nb de frames

            num = nansum( nansum( (Xij - Xi).^2  / (nb_frame*(n_mean-1))));
            den = nansum( nansum( (Xij - X).^2 / (nb_frame*n_mean -1)));

            VR_inter(i) = num / den ;
        end
        All_data(grp).suj_mean(mu).combi_VR_opti = C(min(VR_inter)==VR_inter, :);
        All_data(grp).suj_mean(mu).VR_inter_opti = min(VR_inter);
    end
end

%% kmean sur rmse + corr

suj_set = 1:12;
for grp = 1 : 8
    for mu = 1 : 9
        for k = 1 : 7
            All_data(grp).suj_mean(mu).k(k).idx = kmeans([All_data(grp).suj_mean(mu).rmse All_data(grp).suj_mean(mu).corr] ,k, 'MaxIter',1000, 'Replicates',100);
            All_data(grp).suj_mean(mu).k(k).grp_max = max( unique( sum(All_data(grp).suj_mean(mu).k(k).idx == All_data(grp).suj_mean(mu).k(k).idx')));
        end
        
        n_k = [All_data(grp).suj_mean(mu).k(:).grp_max];
        
        All_data(grp).suj_mean(mu).k_opti = find(min(n_k((n_k-6)>=0)-6) == (n_k((n_k-6)>=0)-6));
        
        All_data(grp).suj_mean(mu).combi_kmeans_opti = suj_set(All_data(grp).suj_mean(mu).k(All_data(grp).suj_mean(mu).k_opti(1)).idx == mode(All_data(grp).suj_mean(mu).k(All_data(grp).suj_mean(mu).k_opti(1)).idx));
    end
end

%%
for grp = 1 : 8
    figure(grp)
    for mu = 1 : 9
        subplot(3,3,mu)
        hold on
        plot(All_data(grp).suj_mean(mu).mu', 'k')
        plot(All_data(grp).suj_mean(mu).mu(All_data(grp).suj_mean(mu).combi_kmeans_opti',:)', 'b') 
        plot(All_data(grp).suj_mean(mu).mu(All_data(grp).suj_mean(mu).combi_VR_opti',:)', 'r')
    end
end
%%

for mu = 1 : 9
    mintersect(All_data(1).suj_mean(mu).combi_VR_opti, All_data(2).suj_mean(mu).combi_VR_opti, All_data(3).suj_mean(mu).combi_VR_opti, All_data(4).suj_mean(mu).combi_VR_opti,...
               All_data(5).suj_mean(mu).combi_VR_opti, All_data(6).suj_mean(mu).combi_VR_opti, All_data(7).suj_mean(mu).combi_VR_opti, All_data(8).suj_mean(mu).combi_VR_opti)
end

for mu = 1 : 9
     a = mintersect(All_data(5).suj_mean(mu).combi_VR_opti, All_data(6).suj_mean(mu).combi_VR_opti) %, ...
%                    All_data(5).suj_mean(mu).combi_VR_opti, All_data(6).suj_mean(mu).combi_VR_opti)
end
