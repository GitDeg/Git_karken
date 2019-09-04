%% Initialisation

clear all;
close all;
clc;

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))

sujet_name={'001','002','003','004','005','006','007','008','009','010','011','012'};

pathname='C:\Users\p1218107\Documents\Data_Piano\data_newREPS'; 
cd(pathname)
%load('All_data.mat')
load('All_data2.mat')


%% trie data chaque sujet 


for grp = 2 % 1 : 8
    for mu = 6 %1 : 9
%         for suj = 1 : 12 
%             %% Variance ratio (coefficient de dispertion) INTRA cycle : entre les cycles d'un meme muscle 
%             % On utilise toutes les données pour calculer la dispertion des cycles
%             % autour de la moyenne. ce coefficient permet de comparer des variations
%             % lorsque les moyennes sont differentes. 
% 
%             nb_frame = length(All_data(grp).suj(suj).mu(mu).data_selec); % nb de frame
%             n_data = size(All_data(grp).suj(suj).mu(mu).data_selec,1) ;
%             
%             Xij = All_data(grp).suj(suj).mu(mu).data_selec;
%             Xi = nanmean(Xij);
%             X = nanmean(nanmean(Xi));
% 
%             num = nansum( nansum( (Xij - Xi).^2  / (nb_frame*(n_data-1))));
%             den = nansum( nansum( (Xij - X).^2 / (nb_frame*n_data -1)));
% 
%             All_data(grp).suj(suj).mu(mu).VR_intra = num / den ;% nb de cycle 
%             
%             
%         end
        
%         for suj = 1 : 12 
%             %% Variance ratio (coefficient de dispertion) INTRA cycle : entre les cycles d'un meme muscle 
%             % On utilise toutes les données pour calculer la dispertion des cycles
%             % autour de la moyenne. ce coefficient permet de comparer des variations
%             % lorsque les moyennes sont differentes. 
% 
%             nb_frame = length(All_data(grp).suj(suj).mu(mu).data_selec_10); % nb de frame
%             n_data = size(All_data(grp).suj(suj).mu(mu).data_selec_10,1) ;
%             
%             Xij = All_data(grp).suj(suj).mu(mu).data_selec_10;
%             Xi = nanmean(Xij);
%             X = nanmean(nanmean(Xi));
% 
%             num = nansum( nansum( (Xij - Xi).^2  / (nb_frame*(n_data-1))));
%             den = nansum( nansum( (Xij - X).^2 / (nb_frame*n_data -1)));
% 
%             All_data(grp).suj(suj).mu(mu).VR_intra_10 = num / den ;% nb de cycle 
%             
%             
%         end
        %% Variance ratio (coefficient de dispertion) INTER individuel : entre les sujets  
        % Cette fois on utilise les données moyennes pour calculer la variation
        % entre les sujet réalisant un meme test. 

        Xij = All_data(grp).suj_mean(mu).mu(:,1:end-1)  ;
        Xi = nanmean(Xij);
        X = nanmean(Xi);

        %n_mean = size(ALL_grp(grp).mu(mu).mean(:,1:id_suj),1);
        n_mean = size(Xij,1); % nb de sujet 
        nb_frame = size(Xij,2); % nb de frames

        num = nansum( nansum( (Xij - Xi).^2))%%  / (nb_frame*(n_mean-1))));
        den = nansum( nansum( (Xij - X).^2))%% / (nb_frame*n_mean -1)));
        
        n_mean = 12
        nb_frame = 629
        
        (nb_frame*n_mean -1)/(nb_frame*(n_mean-1))

        All_data(grp).suj_mean(mu).VR_inter = num / den ;
    end
end


%%
for grp = 1 : 8
    for mu = 1 : 9
        
        nb_frame = length(All_data(grp).suj_mean(mu).mu);
        
        %% Variance ratio (coefficient de dispertion) INTER individuel : entre les sujets  
        % Cette fois on utilise les données moyennes pour calculer la variation
        % entre les sujet réalisant un meme test. 
        for k = 1 : 11
             
            n_suj = [1:12]';
            
            B = unique(All_data(grp).suj_mean(mu).k(k).idx);
            
            Ncount = histc(All_data(grp).suj_mean(mu).k(k).idx, B);
            k_max = find(max(Ncount) == Ncount);
            
            if length(k_max)>1               
                k_max = find(min(All_data(grp).suj_mean(mu).k(k).rmse_k(k_max)) == All_data(grp).suj_mean(mu).k(k).rmse_k);
            end
            
            Xij = All_data(grp).suj_mean(mu).mu(n_suj(All_data(grp).suj_mean(mu).k(k).idx == k_max) ,:)  ;
            Xi = nanmean(Xij);
            X = nanmean(Xi);

            n_mean = size(Xij,1); % nb de sujet 

            num = nansum( nansum( (Xij - Xi).^2  / (nb_frame*(n_mean-1))));
            den = nansum( nansum( (Xij - X).^2 / (nb_frame*n_mean -1)));

            All_data(grp).suj_mean(mu).k(k).VR_inter = num / den ;
            
        end
            
        
    end
end


 %% Essai Variance ratio INTRA_individuel sur fenetre glissante
% 
% for grp = 1 : 8
%     for mu = 1 : 9
%         for suj = 1 : 12
%             
%             k_wind = 50;
%             sujet(suj).grp(grp).mu(mu).VR_intra_gliss = [] ; 
%             
%             for n = 1+k_wind/2 : (id_suj-1) - k_wind/2
%                 % size(sujet(suj).data(mu).grp1(:,1:k),2) ; % nb de frame
%                 n_data = size(sujet(suj).grp(grp).mu(mu).data(:,1:k_wind),1) ;
% 
%                 Xij = sujet(suj).grp(grp).mu(mu).data(:, n-k_wind/2 : n + k_wind/2);
%                 Xi = nanmean(Xij);
%                 X = nanmean(nanmean(Xi));
% 
%                 num = nansum( nansum( (Xij - Xi).^2  / (k_wind*(n_data-1))));
%                 den = nansum( nansum( (Xij - X).^2 / (k_wind*n_data -1)));
% 
%                 sujet(suj).grp(grp).mu(mu).VR_intra_gliss(n) = num / den ;% nb de cycle 
%             end
%             sujet(suj).grp(grp).mu(mu).VR_intra_gliss([1:k_wind/2 n:id_suj-1]) = NaN;
%        end
%     end
% end

%% Variation Coefficient (CV) permit comparison of the variability of data with different mean
% one dimension 

for grp = 1 : 8
    for mu = 1 : 9
        for suj = 1 : 12
            
            sdi = nanstd(All_data(grp).suj(suj).mu(mu).data_selec_10);
            meani = nanmean(All_data(grp).suj(suj).mu(mu).data_selec_10);
            
            All_data(grp).suj(suj).mu(mu).CV = sqrt( mean( sdi.^2)) / mean(abs(meani)) ;
            
            All_data(grp).suj(suj).mu(mu).CVi = sdi./ meani ;
            
        end 
        
        sdi = nanstd(All_data(grp).suj_mean(mu).mu  );
        meani = nanmean(All_data(grp).suj_mean(mu).mu  );
        
        All_data(grp).suj_mean(mu).CV_inter = sqrt( mean( sdi.^2)) / mean(abs(meani)) ;
        
    end
end
%% 

for grp = 1 : 8 
    for mu = 1 : 9 
        for suj = 1 : 12
            GRP(grp).boxplot_data(suj,mu) = All_data(grp).suj(suj).mu(ordre_mu(mu)).CV  ;
            
            MU(mu).boxplot_data(suj,grp) = All_data(grp).suj(suj).mu(ordre_mu(mu)).CV ;
            
            SUJ(suj).boxplot_data(mu,grp) = All_data(grp).suj(suj).mu(ordre_mu(mu)).CV ;
           
        end
%        GRP(grp).boxplot_data_VR_inter(mu) = All_data(grp).suj_mean(ordre_mu(mu)).VR_inter ;
%         MU(mu).boxplot_data_VR_inter_maj(grp) = All_data(grp).suj_mean(mu).VR_inter_maj ;
    end
end
%%
for grp = 1 : 8 
    for mu = 1 : 9
        
        suj_maj = All_data(grp).suj_mean(mu).mu_maj(:,end); 
        
        for suj = 1 : length(suj_maj)
            
            GRP(grp).boxplot_data_maj(suj,mu) = All_data(grp).suj(suj_maj(suj)).mu(mu).VR_intra;
        end
        GRP(grp).boxplot_data_maj(GRP(grp).boxplot_data_maj == 0) = NaN;
    end
end

% 
% subplot(2,1,1);plot(mean(All_data(1).suj(1).mu(4).data_selec_10)); hold on;
% plot(mean(All_data(1).suj(1).mu(4).data_selec_10) + std(All_data(1).suj(1).mu(4).data_selec_10))
% yyaxis right;plot(All_data(1).suj(1).mu(4).CVi')
% 
% subplot(2,1,2);plot(mean(All_data(1).suj(1).mu(9).data_selec_10)); hold on;
% plot(mean(All_data(1).suj(1).mu(9).data_selec_10) + std(All_data(1).suj(1).mu(9).data_selec_10))
% yyaxis right ;plot(All_data(1).suj(1).mu(9).CVi')


% %% muscular synergies for each participant 
% %% Soucis NaN values
% 
% % Que les muscles d'interet / que les sujets aillant tous les muscles 
% % suj = [2:7 9:12]; 
% % mu = [1 2 4 5 6 9] 
% 
% repli = 50 ;
% 
% for syn = 1:5
%     for grp = 1:8
%         for suj = [2:7 9:12]
%         tic
%         data = [sujet(suj).grp(grp).mu([1 2 4 5 6 9]).cycle5]';
%         data(isnan(data(:,1)),:) = [];
%         [sujet(suj).grp(grp).syn(syn).W, sujet(suj).grp(grp).syn(syn).H, sujet(suj).grp(grp).syn(syn).cout] = lee_seung(data,syn,'nbre_iter',repli);
%            
%         sujet(suj).grp(grp).tVAF(syn) = (1- (norm(data - sujet(suj).grp(grp).syn(syn).W* sujet(suj).grp(grp).syn(syn).H,'fro'))^2/...
%                                         (norm(data ,'fro')^2))*100;
%                                     
% %         e = (1 - sum(( data - sujet(suj).grp(grp).syn.W* sujet(suj).grp(grp).syn.H),2).^2./...
% %                                         (sum(data,2).^2))*100;
% %                                     
% %         for mu = 1 : 6         
% %             sujet(suj).grp(grp).VAFm(mu) = e(mu);
% %         end
% 
%         end
%     end
%     save('sujet_Post_ttt.mat','sujet')
% end
% %%
% 

