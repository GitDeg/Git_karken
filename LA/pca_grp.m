%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         Features extractions                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation

clear all;
close all;
clc;

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};

pathname='C:\Users\p1218107\Documents\Data_Piano\data_newREPS\'; 


%% Ouverture Data
% grp = 2; 

Grp_pca = [];
PCA_data = [];
close all
%%

for grp = 5 % : 8
    for suj =  1:12
        Grp_pca(grp).feat(suj).WAMP = zeros(9,630);
        for mu = 1 : 9

            data = All_data(grp).suj_mean(mu).mu(suj,:);
            n_data = length(data);

            data_f = sgolayfilt(data,3,99);

            data_p = diff(data);
            data_p = sgolayfilt(sgolayfilt(data_p,3,99),2,15);
            delta_p = mean(data_p)+std(data_p)*1;

            data_pp = diff(data_p);
            data_pp = sgolayfilt(data_pp,3,99);
            delta_pp = mean(data_pp)+std(data_pp)*2.5;

            %% Mean absolute value (++)

            Grp_pca(grp).feat(suj).MAV(mu) = sum( abs(data))/ n_data;

            %% Waveform length  (++)

            Grp_pca(grp).feat(suj).WL(mu) = sum( abs( data(2:end) - data(1:end-1)));

            %% Willison amplitude (++)

            A = sort(abs(data(1:end-1) - data(2:end)));
            seuil = A(round(n_data*0.95));

            for i = 1 : n_data-1
                if abs(data(i) - data(i+1)) >= seuil
                    f = 1;
                else
                    f = 0;
                end
                Grp_pca(grp).feat(suj).WAMP(mu,i) = f;
            end

             Grp_pca(grp).feat(suj).x_activ(mu,:) = find(Grp_pca(grp).feat(suj).WAMP(mu,:)==1);
    %         eva = evalclusters(x','kmeans','silhouette', 'KList', [1:30]);
    %         feat(suj).n_clust(mu) = sum(diff(feat(suj).x_activ(mu,:))>100) + 1; 
    %         [idx, C] = kmeans(feat(suj).x_activ(mu,:)', feat(suj).n_clust(mu))

    %         figure(mu);
    %         subplot(2,1,1)
    %         plot(data)
    %         subplot(2,1,2)
    %         plot(feat(suj).WAMP(mu,:))
    %    label{mu + (suj-1)*9} = [sujet{suj} '_' All_data(1).suj_mean(mu).name];


            %% Découpage en 3 parties
            x_rep_d = [];
            x_rep_u = [];
            
            [pks,locs] = findpeaks(abs(data_p));
            x_reps = locs(pks>mean(abs(data_p))+std(abs(data_p))*1.2);
            x_reps(x_reps < 125) = [];
            x_reps(x_reps > 500) = [];
            
            for j = 1 : length(x_reps)
                x_rep_d(j) = x_reps(j);
                
                while abs(data_p(x_rep_d(j))) > mean(abs(data_p))/3
                    x_rep_d(j)=x_rep_d(j)-1;
                end
                
                x_rep_u(j) = x_reps(j);
                while abs(data_p(x_rep_u(j))) > mean(abs(data_p))/3
                    x_rep_u(j)=x_rep_u(j)+1;
                end
            end
            
            if isempty(x_rep_d)
                x_rep_d = 1;
            end
            
            if isempty(x_rep_u)
                x_rep_u = length(data);
            end
             Grp_pca(grp).feat(suj).x1(mu) = min(x_rep_d);
             Grp_pca(grp).feat(suj).x2(mu) = max(x_rep_u);
            
%             fig = figure ;
%             plot(data)
            
%             dcm_obj = datacursormode(fig);
%             set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on')
%             
%             B = input('x1 x2  ?');
%             
%             c_info = getCursorInfo(dcm_obj);
%             ind = sort([c_info(:).DataIndex]);
%             
%             Grp_pca(grp).feat(suj).x1(mu) = ind(1);
%             
%             Grp_pca(grp).feat(suj).x2(mu) = ind(2);
            
%             x_rep = NaN(1,length(data_pp));
% 
% 
%             for i = 1:length(data_pp)
%                 if abs(data_p(i)) < delta_p && abs(data_pp(i)) > delta_pp
%                     x_rep(i) = i;
%                 elseif abs(data_p(i)) > delta_p && abs(data_pp(i)) < delta_pp
%                     x_rep(i) = i;
%                 end
%             end
%             x_rep([1:200 end-200:end]) = NaN;
%             x_rep = x_rep(~isnan(x_rep));
%             if isempty(x_rep)
%                 x_rep(1) = 1;
%                 x_rep(2) = length(data_pp);
%             end


            x0 = 1;
%             feat(suj).x1(mu) = x_rep(1);
%             feat(suj).x2(mu) = x_rep(end);
            x3 = length(data);


            %% feature activity initial x0 a x1


            x_deb = x0:Grp_pca(grp).feat(suj).x1(mu);
            y_deb = data(x0:Grp_pca(grp).feat(suj).x1(mu));

            X = [ones(length(x_deb),1) x_deb'];
            Grp_pca(grp).feat(suj).b_deb(mu,:) = X\y_deb';

            y_calc_deb=X*Grp_pca(grp).feat(suj).b_deb(mu,:)';



            %% feature near attack
            
            d_f = sgolayfilt(data(Grp_pca(grp).feat(suj).x1(mu):Grp_pca(grp).feat(suj).x2(mu)),3,15);
 
            [ Grp_pca(grp).feat(suj).mu(mu).pks, Grp_pca(grp).feat(suj).mu(mu).locs, Grp_pca(grp).feat(suj).mu(mu).w, Grp_pca(grp).feat(suj).mu(mu).p] = findpeaks(d_f,'MinPeakProminence',(max(d_f)-min(d_f))/3);
            
            if isempty(Grp_pca(grp).feat(suj).mu(mu).pks)
                Grp_pca(grp).feat(suj).mu(mu).pks = 0;
            end
            if isempty(Grp_pca(grp).feat(suj).mu(mu).locs)
                Grp_pca(grp).feat(suj).mu(mu).locs = 0;
            end
            if isempty(Grp_pca(grp).feat(suj).mu(mu).w)
                Grp_pca(grp).feat(suj).mu(mu).w = 0;
            end
            if isempty(Grp_pca(grp).feat(suj).mu(mu).p)
                Grp_pca(grp).feat(suj).mu(mu).p = 0;
            end
            
            Grp_pca(grp).feat(suj).mu(mu).w_pourc = Grp_pca(grp).feat(suj).mu(mu).w / length(d_f);
    
    
            Grp_pca(grp).feat(suj).pks(mu) =  max(Grp_pca(grp).feat(suj).mu(mu).pks) ;
            Grp_pca(grp).feat(suj).w_pourc(mu) = max(Grp_pca(grp).feat(suj).mu(mu).w_pourc);
            Grp_pca(grp).feat(suj).p(mu) = max(Grp_pca(grp).feat(suj).mu(mu).p);
            Grp_pca(grp).feat(suj).n_pks(mu) = length(Grp_pca(grp).feat(suj).mu(mu).pks);



            %% feature after attack
    %         
    %         feat(suj).p_fin(mu) = (mean(data(x3)) - mean(data(x2)))/(x3 - x2);
            x_fin = Grp_pca(grp).feat(suj).x2(mu):x3;
            y_fin = data(Grp_pca(grp).feat(suj).x2(mu):x3);

            X = [ones(length(x_fin),1) x_fin'];
            Grp_pca(grp).feat(suj).b_fin(mu,:) = X\y_fin';

            y_calc_fin=X*Grp_pca(grp).feat(suj).b_fin(mu,:)';

    %%
            figure(mu)
            subplot(4,3,suj)
            plot(data)
            hold on 
            line([Grp_pca(grp).feat(suj).x1(mu) Grp_pca(grp).feat(suj).x1(mu)], ylim)
            line([Grp_pca(grp).feat(suj).x2(mu) Grp_pca(grp).feat(suj).x2(mu)], ylim)
            plot(x_fin, y_calc_fin, 'r')
            plot(x_deb, y_calc_deb, 'r')

        end

    %    PCA_data = [PCA_data; feat(suj).MAV feat(suj).WL feat(suj).n_clust]; 
    %     PCA_data = [PCA_data; feat(suj).p_ini feat(suj).std_att feat(suj).tau_acti feat(suj).p_fin]; 
    %     Z = zscore (PCA_data);
    %     mapcaplot(Z,sujet)
    end
end

%%

grp = 5;
for mu = 7 %1 : 9
    PCA_data = [];
    for suj = 1:12
        %% PCA
        PCA_data = [PCA_data; Grp_pca(grp).feat(suj).x1(mu) Grp_pca(grp).feat(suj).x2(mu) Grp_pca(grp).feat(suj).b_deb(mu,:) Grp_pca(grp).feat(suj).b_fin(mu,:) ...
            Grp_pca(grp).feat(suj).pks(mu) Grp_pca(grp).feat(suj).w_pourc(mu) Grp_pca(grp).feat(suj).p(mu) Grp_pca(grp).feat(suj).n_pks(mu)]; 

    end
    Z = zscore (PCA_data);
    mapcaplot(Z,sujet)
end

%%

idx = kmeans(Z,4);

%%

grp = 5;

close all

for mu = 1 : 9
    figure(mu)
    
    for suj = 1 : 12
        subplot(3,4,suj)
        plot(All_data(grp).suj_mean(mu).mu(suj,:))
    end
    suptitle(All_data(1).suj_mean(mu).name)
    
end


%%
grp = 1;
suj = 1; 
mu = 8;
data = All_data(grp).suj_mean(mu).mu(suj,:);

figure; 
plot(data)
hold on 
line(xlim, [mean(data) mean(data)])

line([250 250],ylim)
line([400 400],ylim)


%% 


col = {'g' 'b' 'y' 'm' 'c' 'k'};
suj = 1 : 12 ; 

for grp = 1 : 8
    close all
    for mu = 1 : 9

        A = corrcoef(All_data(grp).suj_mean(mu).mu') ;

%         Z = zscore (All_data(grp).suj_mean(mu).rmse);
        Z = zscore ([All_data(grp).suj_mean(mu).rmse A]);
        %mapcaplot(Z)

%        [coeff,score,latent] = pca(Z);
        
        %% Retrait des sujets aberrants 
        
        Ncount = 1 ; 
        
        while min(Ncount) == 1
            
            OptiK = 2;        
%             clust_opti = kmeans(score(:,1:3),OptiK(1), 'emptyaction', 'singleton', 'replicate', 50);
            clust_opti = kmeans(Z,OptiK(1), 'emptyaction', 'singleton', 'replicate', 50);

    %         Ncount_diff = find(diff(Ncount_max)<-1);
    %         %Ncount_diff(Ncount_diff == 1) = [];
    %         if isempty(Ncount_diff)
    %             Ncount_diff = 7;
    %         end
    %         va = evalclusters(score(:,1:3), clust, 'CalinskiHarabasz');

    %         OptiK = Ncount_diff(1) + 1;
            Ncount = hist(clust_opti(~isnan(clust_opti)), unique(clust_opti(~isnan(clust_opti))));


            N_count_max = (find( max(Ncount) == Ncount));
            N_count_min = (find( min(Ncount) == Ncount));
            
            if min(Ncount) == 1
                
%                 score(clust_opti == N_count_min,:) = NaN ;  
                Z(clust_opti == N_count_min,:) = NaN ; 
                
            end
        end
        
        All_data(grp).suj_mean(mu).mu_ssAberrant = [All_data(grp).suj_mean(mu).mu( ~isnan(clust_opti),:) suj(~isnan(clust_opti))'];
        %% Determination K opti : groupe max d'au moins 7 sujets
% 
%         clust = [];
%         Ncount_max = []; 
%         Ncount_diff = []; 
%         
%         for i = 1 : 5
%             clust(:,i) = kmeans(score(:,1:3),i, 'emptyaction', 'singleton', 'replicate', 50);
%             Ncount_max = [Ncount_max max(hist(clust(:,i), unique(clust(:,i))))];
%         end
%         
%         OptiK = find(min(Ncount_max(Ncount_max >= 7)) == Ncount_max);
%         OptiK = OptiK(end);
%         
%         clust_opti = kmeans(score(:,1:3),OptiK(1), 'emptyaction', 'singleton', 'replicate', 50);
        
%         Ncount_diff = find(diff(Ncount_max)<-1);
%         %Ncount_diff(Ncount_diff == 1) = [];
%         if isempty(Ncount_diff)
%             Ncount_diff = 7;
%         end
%         va = evalclusters(score(:,1:3), clust, 'CalinskiHarabasz');

%         OptiK = Ncount_diff(1) + 1;
%         Ncount = hist(clust_opti, unique(clust_opti));
%         N_count_max = (find( max(Ncount) == Ncount));
        
        

%% Determination du k opti : min VR_inter

%         OptiK = find( min([All_data(grp).suj_mean(mu).k(1:5).VR_inter]) == [All_data(grp).suj_mean(mu).k(1:5).VR_inter]);
        


%         if length(OptiK)>1
%             OptiK = OptiK(find( min([All_data(grp).suj_mean(mu).k(OptiK).rmse_mean]) == [All_data(grp).suj_mean(mu).k(OptiK).rmse_mean]));
%         end
%         
%         clust_opti = kmeans(score(:,1:3),OptiK(1), 'emptyaction', 'singleton', 'replicate', 50);
%         Ncount = hist(clust_opti, unique(clust_opti));
%         
%         N_count_max = (find( max(Ncount) == Ncount));

%% Plot 
        
        figure(mu+9*(grp-1))
        subplot(2,2,1)
        plot(All_data(grp).suj_mean(mu).mu')
        xlim([0 length(All_data(grp).suj_mean(mu).mu)])
        
        subplot(2,2,3)
        plot(All_data(grp).suj_mean(mu).mu_ssAberrant') %, 'color', [0 0 0 0.5])
        hold on
%         plot(All_data(grp).suj_mean(mu).mu( clust_opti == N_count_max(1),:)','LineWidth',1.5,'Color', 'r')
        xlim([0 length(All_data(grp).suj_mean(mu).mu)])
        legend(num2str(sum(Ncount)))

        subplot(1,2,2)
        plot(All_data(grp).suj_mean(mu).mu', 'color', [0 0 0 0.5])
        hold on
        plot(All_data(grp).suj_mean(mu).mu(isnan(clust_opti),:)','LineWidth',1.5,'Color', 'r')

        legend(num2str(12 - sum(Ncount)))
   
        set(gcf, 'Position', get(0, 'Screensize'))
        suptitle(['grp' num2str(grp) ' / ' All_data(grp).suj_mean(mu).name])
        
        saveas(figure(mu+9*(grp-1)), ['C:\Users\p1218107\Documents\Data_Piano\graphiques\New_reps\Clusters_kmean_pca\','grp' num2str(grp) '_' All_data(grp).suj_mean(mu).name '_ssAberrant2'],'jpg')
    
%% Save des groupes maj

%         All_data(grp).suj_mean(mu).mu_ssAberrant= [All_data(grp).suj_mean(mu).mu( clust_opti == N_count_max(1),:) suj(clust_opti == N_count_max(1))'] ;
    
    end
    
end
%%
figure;
scatter3(score(:,1),score(:,2),score(:,3),100*ones(12,1), clust_opti(:,va.OptimalK))
axis equal


%% Auto-regressive coeffivients (++)

% AR

%% Mean absolute value slope (++)

% MAVS


%%

% %% Integrated EMG
% 
% IEMG = sum( abs(data));
% 
% 
% %% Modified mean absolute value type I
% 
% s_MAV1 = 0;
% 
% for i = 1 : n_data
%     if i <= 0.75*n_data && i >= 0.25 * n_data
%         w = 1;
%     else
%         w = 0.5;
%     end
%     s_MAV1 = s_MAV1 + w*abs( data(i));
% end
% 
% MAV1 = s_MAV1/n_data;
% 
% %% Modified mean absolute value type 2
% 
% s_MAV2 = 0;
% 
% for i = 1 : n_data
%     if i <= 0.75*n_data && i >= 0.25 * n_data
%         w = 1;
%     elseif i < 0.25*n_data
%         w = 4*i/n_data;
%     else
%         w = 4*(i-n_data)/n_data;
%     end
%     s_MAV2 = s_MAV2 + w*abs( data(i));
% end
% 
% MAV2 = s_MAV2/ n_data;
% 
% %% Simple square integral
% 
% SSI = sum( data.^2);
% 
% %% Variance of EMG
% 
% VAR = std(data)^2;
% 
% %% Absolute value of the 3rd, 4rth and 5th temporal moment
% 
% TM3 = abs( sum(data.^3)/n_data);
% 
% TM4 = sum(data.^4)/n_data;
% 
% TM5 = abs( sum(data.^5)/n_data);
% 
% %% Root mean square 
% 
% RMS = sqrt( sum( data.^2)/n_data);
% 
% %% Log detector 
% 
% LOG = exp(sum(log(abs(data)))/n_data);
% 
% 
% %% Average amplitude change 
% 
% AAC = sum( abs( data(2:end) - data(1:end-1))) / n_data;
% 
% %% Difference absolute standard deviation value
% 
% DASDV = sqrt( sum( (data(2:end) - data(1:end-1)).^2 )/(n_data-1));


