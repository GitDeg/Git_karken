%% Initialisation

clear all;
close all;
clc;

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))

sujet_name={'001','002','003','004','005','006','007','008','009','010','011','012'};

GrpNames={'Membre supérieur pressé staccato',...
              'Membre supérieur frappé staccato',...
              'Bassin pressé staccato',...
              'Bassin frappé staccato',...
              'Membre supérieur pressé tenue',...
              'Membre supérieur frappé tenue',...
              'Bassin pressé tenue'...
              'Bassin frappé tenue'};

pathname='C:\Users\p1218107\Documents\Data_Piano\data_newREPS'; 
cd(pathname)
load('All_data.mat')



%%

% sujet = [];
% sujet.name = [];
% sujet.grp = [];
% 
% 
% sujet.grp.name = [];
% sujet.grp.mu = [] ; 
% 
% sujet.grp.mu.name = []; 
% sujet.grp.mu.data = [];
% sujet.grp.mu.mean = [];
% sujet.grp.mu.cycle5 = []; 
% sujet.grp.mu.VR_intra = [];
% sujet.grp.mu.VR_inter = [];
% 
% sujet.grp.syn = [];
% sujet.grp.syn.W = [];
% sujet.grp.syn.H = [];
% sujet.grp.syn.cout = [];
% 
% sujet.grp.tVAF = [];
% sujet.grp.VAFm = [];

%% trie data chaque sujet 


for grp = 1 : 8
    for mu = 1 : 10
        for suj = 1 : 12
            
            sujet(suj).name = sujet_name{suj};
            
            sujet(suj).grp(grp).name = GrpNames{grp};
            
            %sujet(suj).grp(grp).mu(mu).name = ALL_grp(1).mu(mu).name; 
            sujet(suj).grp(grp).mu(mu).name = Post_ttt(1).suj_mean(mu).name;
            
            %sujet(suj).grp(grp).mu(mu).data = ALL_grp(grp).mu(mu).data(ALL_grp(grp).mu(mu).data(:,id_suj)==suj,:) ; 
            sujet(suj).grp(grp).mu(mu).data = Post_ttt(grp).suj(suj).mu(mu).data ;  
            
            %sujet(suj).grp(grp).mu(mu).mean = ALL_grp(grp).mu(mu).mean(ALL_grp(grp).mu(mu).mean(:,id_suj)==suj,:) ;
            sujet(suj).grp(grp).mu(mu).mean = mean(Post_ttt(grp).suj(suj).mu(mu).data);
            %% 5 cycles bout a bout
            
            sujet(suj).grp(grp).mu(mu).data(isnan(sujet(suj).grp(grp).mu(mu).data(:,1)), :)=[];
            if size(sujet(suj).grp(grp).mu(mu).data,1)>nb_cycle
                ran1 = randperm(size(sujet(suj).grp(grp).mu(mu).data,1));
                sujet(suj).grp(grp).mu(mu).cycle5 = reshape(sujet(suj).grp(grp).mu(mu).data(ran1(1:nb_cycle),1:id_suj-1)',[],1);
            else
                sujet(suj).grp(grp).mu(mu).cycle5 = NaN(nb_cycle*(id_suj-1),1);
            end
            %% Variance ratio (coefficient de dispertion) INTRA cycle : entre les cycles d'un meme muscle 
            % On utilise toutes les données pour calculer la dispertion des cycles
            % autour de la moyenne. ce coefficient permet de comparer des variations
            % lorsque les moyennes sont differentes. 

            k = id_suj-1; % size(sujet(suj).data(mu).grp1(:,1:k),2) ; % nb de frame
            n_data = size(sujet(suj).grp(grp).mu(mu).data(:,1:k),1) ;
            
            Xij = sujet(suj).grp(grp).mu(mu).data(:,1:k);
            Xi = nanmean(Xij);
            X = nanmean(nanmean(Xi));

            num = nansum( nansum( (Xij - Xi).^2  / (k*(n_data-1))));
            den = nansum( nansum( (Xij - X).^2 / (k*n_data -1)));

            sujet(suj).grp(grp).mu(mu).VR_intra = num / den ;% nb de cycle 
            
            %% Variance ratio (coefficient de dispertion) INTER individuel : entre les sujets  
            % Cette fois on utilise les données moyennes pour calculer la variation
            % entre les sujet réalisant un meme test. 
            
            %Xij = ALL_grp(grp).mu(mu).mean(:,1:id_suj);
            Xij = Post_ttt(grp).suj_mean(mu).data;
            Xi = nanmean(Xij);
            X = nanmean(Xi);
   
            k = id_suj-1 ; % nb de frame
            %n_mean = size(ALL_grp(grp).mu(mu).mean(:,1:id_suj),1);
            n_mean = size(Post_ttt(1).suj_mean(1).data,1); % nb de sujet 

            num = nansum( nansum( (Xij - Xi).^2  / (k*(n_mean-1))));
            den = nansum( nansum( (Xij - X).^2 / (k*n_mean -1)));
    
            sujet(suj).grp(grp).mu(mu).VR_inter = num / den ;
            
        end
    end
end

%% Essai Variance ratio INTRA_individuel sur fenetre glissante

for grp = 1 : 8
    for mu = 1 : 10
        for suj = 1 : 12
            
            k_wind = 50;
            sujet(suj).grp(grp).mu(mu).VR_intra_gliss = [] ; 
            
            for n = 1+k_wind/2 : (id_suj-1) - k_wind/2
                % size(sujet(suj).data(mu).grp1(:,1:k),2) ; % nb de frame
                n_data = size(sujet(suj).grp(grp).mu(mu).data(:,1:k_wind),1) ;

                Xij = sujet(suj).grp(grp).mu(mu).data(:, n-k_wind/2 : n + k_wind/2);
                Xi = nanmean(Xij);
                X = nanmean(nanmean(Xi));

                num = nansum( nansum( (Xij - Xi).^2  / (k_wind*(n_data-1))));
                den = nansum( nansum( (Xij - X).^2 / (k_wind*n_data -1)));

                sujet(suj).grp(grp).mu(mu).VR_intra_gliss(n) = num / den ;% nb de cycle 
            end
            sujet(suj).grp(grp).mu(mu).VR_intra_gliss([1:k_wind/2 n:id_suj-1]) = NaN;
       end
    end
end

%% Variation Coefficient (CV) permit comparison of the variability of data with different mean
% one dimension 

for grp = 1 : 8
    for mu = 1 : 10
        for suj = 1 : 12
            
            sdi = nanstd(sujet(suj).grp(grp).mu(mu).data(:,1:end-1));
            meani = nanmean(sujet(suj).grp(grp).mu(mu).data(:,1:end-1));
            
            sujet(suj).grp(grp).mu(mu).CV = sqrt( mean( sdi.^2)) / mean(meani) ;
            
            sujet(suj).grp(grp).mu(mu).CVi = sdi./ meani ;
            
        end 
    end
end

%% muscular synergies for each participant 
%% Soucis NaN values

% Que les muscles d'interet / que les sujets aillant tous les muscles 
% suj = [2:7 9:12]; 
% mu = [1 2 4 5 6 9] 

repli = 50 ;

for syn = 1:5
    for grp = 1:8
        for suj = [2:7 9:12]
        tic
        data = [sujet(suj).grp(grp).mu([1 2 4 5 6 9]).cycle5]';
        data(isnan(data(:,1)),:) = [];
        [sujet(suj).grp(grp).syn(syn).W, sujet(suj).grp(grp).syn(syn).H, sujet(suj).grp(grp).syn(syn).cout] = lee_seung(data,syn,'nbre_iter',repli);
           
        sujet(suj).grp(grp).tVAF(syn) = (1- (norm(data - sujet(suj).grp(grp).syn(syn).W* sujet(suj).grp(grp).syn(syn).H,'fro'))^2/...
                                        (norm(data ,'fro')^2))*100;
                                    
%         e = (1 - sum(( data - sujet(suj).grp(grp).syn.W* sujet(suj).grp(grp).syn.H),2).^2./...
%                                         (sum(data,2).^2))*100;
%                                     
%         for mu = 1 : 6         
%             sujet(suj).grp(grp).VAFm(mu) = e(mu);
%         end

        end
    end
    save('sujet_Post_ttt.mat','sujet')
end
%%


