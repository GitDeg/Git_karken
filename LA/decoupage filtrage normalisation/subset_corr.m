k_min = [];

for grp = 1% : 8
    for mu = 1 : 9

        k_min = [k_min, size(All_data(grp).suj(suj).mu(mu).data_norm,1)];

    end
end



%%
for suj = 1 : 12
    for grp = 1 : 8
        for mu = 1 : 9
            tic
            All_data(grp).suj(suj).mu(mu).data_selec = [];
            DATA = [mean(All_data(grp).suj(suj).mu(mu).data_norm)', All_data(grp).suj(suj).mu(mu).data_norm']; 
            
            R = corrcoef(DATA);

            R_tot = sort(R(1,:));
            
            n = [];
            
            for k = 1 : 7
                n(k) = find( R_tot(end - k)== R(1,:));
            end
            
            All_data(grp).suj(suj).mu(mu).data_selec = [DATA(:,  n)', ones(7,1)*suj];
%             v = 1 : size(DATA,2);
%             k = 5;
% 
%             C = nchoosek(v,k);
%             R_tot = zeros(length(C),10);
%             
%             for i = 1 : length(C)
%                 C_tot =  nchoosek(C(i,:),2);
%                  
%                 for j = 1 : length(C_tot)
%                     R_tot(i,j)= R(C_tot(j,1),C_tot(j,2));
%                 end
%                 
%             end
%             toc
%             All_data(grp).suj(suj).mu(mu).data_selec = DATA(:,  C( max( mean( R_tot,2)) == mean( R_tot,2), :))';
            
        end
    end
end

%%

figure ;
subplot(2,1,1)
plot(All_data(grp).suj(suj).mu(mu).data_norm') ;

subplot(2,1,2)
plot(All_data(grp).suj(suj).mu(mu).data_selec');

%% 

for grp = 1 : 8
    for mu = 1 : 9
        
        All_data(grp).suj_mean(mu).mean_selec = NaN(12,630);
        DATA = [mean(All_data(grp).suj_mean(mu).mu)', All_data(grp).suj_mean(mu).mu']; 

        R = corrcoef(DATA);

        R_tot = sort(R(1,:));

        n = [];

        for k = 1 : 10
            n(k) = find( R_tot(end - k)== R(1,:));
        end

        All_data(grp).suj_mean(mu).mean_selec(n-1,:) = DATA(:,  n)';
%             v = 1 : size(DATA,2);
%             k = 5;
% 
%             C = nchoosek(v,k);
%             R_tot = zeros(length(C),10);
%             
%             for i = 1 : length(C)
%                 C_tot =  nchoosek(C(i,:),2);
%                  
%                 for j = 1 : length(C_tot)
%                     R_tot(i,j)= R(C_tot(j,1),C_tot(j,2));
%                 end
%                 
%             end
%             toc
%             All_data(grp).suj(suj).mu(mu).data_selec = DATA(:,  C( max( mean( R_tot,2)) == mean( R_tot,2), :))';

    end
end

%% Correlation

for grp = 1 : 8
    for mu = 1 : 9
        
            tic

            DATA = All_data(grp).suj_mean(mu).mu'; 
            R = corrcoef(DATA);

            v = 1 : size(DATA,2);
            
        for k = 2:12
            
            n = [];
            C = nchoosek(v,k);
            R_tot = [];

%             for i = 1 : size(C,1)
%                 C_tot =  nchoosek(C(i,:),2);
% 
%                 for j = 1 : size(C_tot,1)
%                     R_tot(i,j)= R(C_tot(j,1),C_tot(j,2));
%                 end
% 
%             end
%             toc
%             corr_grp(grp).mu(mu).k(k).corr_max = R_tot( find( max( mean( R_tot,2)) == mean( R_tot,2)),:);
%             
%             corr_grp(grp).mu(mu).k(k).mean_corr_max = max( mean( R_tot,2)) ;
%             
%             corr_grp(grp).mu(mu).k(k).combi_max = C(find( max( mean( R_tot,2)) == mean( R_tot,2)),:);
            
            
            
            
            
        end
        
    end
end

%% Plot correlation

for grp = 1 : 8
    for mu = 1 : 9
        figure(mu)
        for i = 2 : length(corr_grp(grp).mu(mu).k  )
            subplot(4,3,i)
            plot(All_data(grp).suj_mean(mu).mu( corr_grp(grp).mu(mu).k(i).combi_max',:)')
            legend(num2str( corr_grp(grp).mu(mu).k(i).mean_corr_max))
        end
    end
end
    
%% RMSE
for grp = 1 : 8
    for mu = 1 : 9
            DATA = All_data(grp).suj_mean(mu).mu'; 
            
        for i = 1 : size(DATA,2)
            for j = 1 : size(DATA,2)
                All_data(grp).suj_mean(mu).rmse(i,j) = sqrt( sum( (DATA(:,i) - DATA(:,j)).^2))/numel(DATA(:,i));
            end
        end
        
    end
end
%% Plot correlation 


for grp = 1 : 8
    for mu = 1 : 9
        for k = 1 : 12
            All_data(grp).suj_mean(mu).k(k).idx = kmeans(All_data(grp).suj_mean(mu).rmse,k);
        end
    end
end

%%

n_suj = 1 : 12;
for grp = 1 : 8
    for mu = 1 : 9
        for k = 1 : 12
            for i = 1 : k

                n_suj_k = n_suj(All_data(grp).suj_mean(mu).k(k).idx == i);

                if length(n_suj_k) > 1

                    C_tot = nchoosek(n_suj_k,2);
                    rmse_sum = [];

                    for j = 1 : size(C_tot,1)
                        rmse_sum(j) = All_data(grp).suj_mean(mu).rmse(C_tot(j,1),C_tot(j,2));
                    end
                    All_data(grp).suj_mean(mu).k(k).rmse_k(i) = mean(rmse_sum);

                else
                    All_data(grp).suj_mean(mu).k(k).rmse_k(i) = 0;
                end
            end
            All_data(grp).suj_mean(mu).k(k).rmse_mean = mean(All_data(grp).suj_mean(mu).k(k).rmse_k);
        end
    end
end

%%
figure;
plot(DATA(:,idx==1),'b')

hold on

plot(DATA(:,idx==2),'r')

plot(DATA(:,idx==3),'g')
plot(DATA(:,idx==4),'m')



