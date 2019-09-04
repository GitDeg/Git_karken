%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     Normalisation et filtrage stats                     % 
%                       par moyenne des conditions                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation

clear all;
close all;
clc;

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};

pathname='C:\Users\p1218107\Documents\Data_Piano\data_newREPS\'; 

alph = 0.01;


%% Ouverture des fichiers 

for suj = 10 %1 : 12
    
    pathname_emg = [pathname,sujet{suj},'_data.mat'];
    load(pathname_emg)
    
    for grp = 1 : 8
        for mu = 1 : 9
            
%        Test(grp).suj(suj).mu(mu).muscle = GRP(grp).cycle(mu).muscle;
        Test(grp).suj(suj).mu(mu).data = GRP(grp).cycle(mu).data; 
        
        end 
    end
end


%%
%% criteria III 
% 50% different from others and 20% outlying 
% Test(grp).suj(suj).mu(mu).out_pourc = 100*sum(Test(grp).suj(suj).mu(mu).out,2)/sum(sum(Test(grp).suj(suj).mu(mu).out,2));

for grp = 1 : 8
    for suj = 1 : 12
        for mu = 1 : 9 
            
           
            Test(grp).suj(suj).mu(mu).data_ttt = Test(grp).suj(suj).mu(mu).data;
            
            while 1
                Test(grp).suj(suj).mu(mu).outliers = [];
                Test(grp).suj(suj).mu(mu).diff = [];
                Test(grp).suj(suj).mu(mu).diff_pourc  = [];
                Test(grp).suj(suj).mu(mu).normality = [];
                Test(grp).suj(suj).mu(mu).out = [];
                Test(grp).suj(suj).mu(mu).out = [];
                
                Test(grp).suj(suj).mu(mu).outliers = zeros(1,size(Test(grp).suj(suj).mu(mu).data_ttt,1));
                
                Y = Test(grp).suj(suj).mu(mu).data_ttt;
                
                n1 = size(Y,1);
            
                [p, tbl, stats] = friedman(Y',1, 'off');
                c = multcompare(stats, 'Alpha', alph, 'Display', 'off'); 

                for j = 1 : size(Y,1)
                    dif = c( logical( nansum( c(:,1:2) == j, 2)), :);
                    Test(grp).suj(suj).mu(mu).diff(j) = length( dif( dif(:,6) < alph)) ;
                end
                Test(grp).suj(suj).mu(mu).diff_pourc = 100*(Test(grp).suj(suj).mu(mu).diff)/size(Test(grp).suj(suj).mu(mu).data,1);

                for i = 1 : length(Y)
                   Test(grp).suj(suj).mu(mu).normality(i) = kstest( Y( :,i));
                   Test(grp).suj(suj).mu(mu).out(:,i) = isoutlier( Y( :,i)); 
                end
                Test(grp).suj(suj).mu(mu).out_pourc = 100*sum(Test(grp).suj(suj).mu(mu).out,2)/sum(sum(Test(grp).suj(suj).mu(mu).out,2));
                
                for cycl = 1 : length(Test(grp).suj(suj).mu(mu).diff_pourc) 
                
%                     if Test(grp).suj(suj).mu(mu).diff_pourc(cycl) > 30 &&  Test(grp).suj(suj).mu(mu).out_pourc(cycl) > 30
%                         Test(grp).suj(suj).mu(mu).outliers(cycl) = 1;
%                     end
% 
%                     if Test(grp).suj(suj).mu(mu).diff_pourc(cycl) > 40 &&  Test(grp).suj(suj).mu(mu).out_pourc(cycl) > 25
%                         Test(grp).suj(suj).mu(mu).outliers(cycl) = 1;
%                     end

                    if Test(grp).suj(suj).mu(mu).diff_pourc(cycl) > 50 &&  Test(grp).suj(suj).mu(mu).out_pourc(cycl) > 20
                        Test(grp).suj(suj).mu(mu).outliers(cycl) = 1;
                    end

                    if Test(grp).suj(suj).mu(mu).diff_pourc(cycl) > 60 &&  Test(grp).suj(suj).mu(mu).out_pourc(cycl) > 15
                        Test(grp).suj(suj).mu(mu).outliers(cycl) = 1;
                    end
% 
%                     if Test(grp).suj(suj).mu(mu).diff_pourc(cycl) > 70 &&  Test(grp).suj(suj).mu(mu).out_pourc(cycl) > 10
%                         Test(grp).suj(suj).mu(mu).outliers(cycl) = 1;
%                     end
% 
%                     if Test(grp).suj(suj).mu(mu).diff_pourc(cycl) > 80 &&  Test(grp).suj(suj).mu(mu).out_pourc(cycl) > 5
%                         Test(grp).suj(suj).mu(mu).outliers(cycl) = 1;
%                     end

                end   
                
                Test(grp).suj(suj).mu(mu).data_ttt( logical( Test(grp).suj(suj).mu(mu).outliers'),:) = [];
                
                n2 = size(Test(grp).suj(suj).mu(mu).data_ttt,1);
                %n2 = size(Y,1) - sum( double(isnan( Test(grp).suj(suj).mu(mu).data_ttt(:,1))));
                
                if n2 == n1
                    break
                end
            end
%             Test(grp).suj(suj).mu(mu).outliers_tot = sum(Test(grp).suj(suj).mu(mu).outliers);
        end
    end
end

%% Retrait manuel des cycles les plus improbables 
% ils peuvent grandement influencer la normalisation

for grp = 1 : 8
    for suj = 10% 1 : 12
        suj
        for mu = 1 : 9 
            
             Save = 'N';
             A = Test(grp).suj(suj).mu(mu).data; %_ttt ;
             while Save == 'N'
                close all 
             
                 fig = figure ;
                 suptitle([sujet{suj} ' mu ' num2str(mu)])
                 subplot(2,1,1)
                 plot(A')
                 hold on
                 plot(mean(A ) + 3*std(A),'linewidth',2 ,'color', 'r')
                 n_out = [] ; 
                 
                 modif = input('Modifications ? [ Y / N ]','s');
                 if modif ~= 'N'
                     dcm_obj = datacursormode(fig);
                     set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on')
                     
                     str='N';
                     while str == 'N'
                            str = input('Surligner les courbes ? [ Y / N ]','s');
                            
                            if str ~= 'N'
                                c_info = getCursorInfo(dcm_obj);
                                hold on
                                for i = 1: length(c_info)
                                    n_out(i) = find(A(:,c_info(i).DataIndex) == max(A(:,c_info(i).DataIndex)));
                                    hold on
                                    plot(A(n_out(i),:)', 'linewidth', 1.5, 'color', 'g')
                                end
                                
                                A(n_out,:) = [];
                                subplot(2,1,2)
                                plot(A')
                                hold on
                                plot(mean(A ) + 3*std(A),'linewidth',2 ,'color', 'r')
                            end
                     end
                 end
                 
                 Save = input('Save ? [ Y / N ] ', 's');
                 if Save ~= 'N'
                     %Test(grp).suj(suj).mu(mu).data_ttt(n_out,:) = [];
                     Test(grp).suj(suj).mu(mu).data_ttt = A;
                 end
                 
             end
            
        end
    end
end



%%
for grp = 1 : 8
    for suj = 1 : 12
        figure(12*(grp-1) + suj)
        for mu = 1 : 9   
%             if size(Test(grp).suj(suj).mu(mu).data,1) ~= size(Test(grp).suj(suj).mu(mu).data_ttt,1)
                %figure(suj) 
                %suptitle([Test(grp).suj(suj).mu(mu).muscle ' ' sujet{suj} ' ' GRP(grp).name 'dif =' num2str(abs(size(Test(grp).suj(suj).mu(mu).data,1) - size(Test(grp).suj(suj).mu(mu).data_ttt,1)))] )
%                 subplot(3,3,mu)
%                 plot(Test(grp).suj(suj).mu(mu).data')
                suptitle([ 'suj = ' sujet{suj} ' grp = ' num2str(grp)])
                subplot(3,3,mu)
                plot(Test(grp).suj(suj).mu(mu).data_ttt')
%             end
        end
    end
end

%%
% suj = 11;
% grp = 8;
% mu = 5;
% 
% figure
% subplot(2,1,1)
% plot(Test(grp).suj(suj).mu(mu).data_ttt')
% 
% x_out = 1775; 
% 
% n = find(Test(grp).suj(suj).mu(mu).data_ttt(:,x_out) == max(Test(grp).suj(suj).mu(mu).data_ttt(:,x_out)));
% 
% hold on
% plot(Test(grp).suj(suj).mu(mu).data_ttt(n,:)', 'linewidth', 2, 'color', 'r')
% 
% Test(grp).suj(suj).mu(mu).data_ttt(n,:) = [];
% 
% subplot(2,1,2)
% plot(Test(grp).suj(suj).mu(mu).data_ttt')
for grp = 1 : 8
    All_data(grp).name = Test(grp).name;
end
%%

% pathname_emg = [pathname,'test_ttt.mat'];
% load(pathname_emg)

for suj = 9 %1 : length(sujet)
    
    for mu = 4 %1 : 9 %length(Test(1).suj(1).mu)
        EMG_mu=[];
        
        for grp = 1 : 8 %length(Test)
            EMG_mu=[EMG_mu; Test(grp).suj(suj).mu(mu).data_ttt];
        end
        
        EMG_mean=nanmean(nanmean(EMG_mu));
        
        for grp = 1 : 8  
            All_data(grp).suj_mean(mu).name = Test(grp).suj(suj).mu(mu).muscle;
            All_data(grp).suj_mean(mu).mu(suj,:) = downsample( mean( Test(grp).suj(suj).mu(mu).data_ttt / EMG_mean), 5);
            All_data(grp).suj(suj).mu(mu).name = Test(grp).suj(suj).mu(mu).muscle;
            All_data(grp).suj(suj).mu(mu).data_norm =  downsample(Test(grp).suj(suj).mu(mu).data_ttt' / EMG_mean,5)';
        end
    end 
end

%%

for suj = 1 : 12
    for mu = 1 : 9 %length(Test(1).suj(1).mu)
        for grp = 1 : 8  
            
            All_data(grp).suj_mean(mu).mu(suj,:) = mean( All_data(grp).suj(suj).mu(mu).data_selec_10);

        end
    end 
end
%%

close all

for mu = 1 : 9
    figure(mu)
    for grp = 1 : 8
        subplot(4,2,grp)
        plot(All_data(grp).suj_mean(mu).mu')
        hold on 
%        plot(All_data(grp).suj_mean(mu).mu(10,:), 'r', 'linewidth', 2)
%        plot(All_data(grp).suj_mean(mu).mean_selec')
    end
    suptitle(All_data(grp).suj_mean(mu).name)
end


%%
% 
% for mu = 1 : 10 
%     for grp = 1 : 8 
%         
% %         A = [];
% %         
% %         
% %         for suj = 1 : 12
% %             A = [A; Test(grp).suj(suj).mu(mu).mean]; 
% %         end
% %         
%         A = Test(grp).suj_mean(mu).data ; 
%         
%         Xij = Test(grp).suj_mean(mu).data(:,1:end);
%         Xi = nanmean(Xij);
%         X = nanmean(Xi);
% 
%         n_mean = size(Test(grp).suj_mean(mu).data,1); % nb de sujet 
% 
%         num = nansum( nansum( (Xij - Xi).^2  / (788*(n_mean-1))));
%         den = nansum( nansum( (Xij - X).^2 / (788*n_mean -1)));
% 
%         Test(grp).suj_mean(mu).VR_inter = num / den ;
%             
%         figure(mu)
%         subplot(4,2,grp)
%         spm1d.plot.plot_meanSD(1:length(A),A)
%         hold on 
%         plot(A')
%         legend(num2str(Test(grp).suj_mean(mu).VR_inter))
%         
%     end
% end
% %%
% B = [];
% close all
% 
% for mu = 1 : 10
%     for grp = 1 : 8 
%   
%         Test(grp).suj_mean(mu).outliers = zeros(1,size(Test(grp).suj_mean(mu).data,1));
%         Test(grp).suj_mean(mu).data_ttt = Test(grp).suj_mean(mu).data;
% 
% 
%             
%         if Test(grp).suj_mean(mu).VR_inter > 0.75
%         
%             Y = Test(grp).suj_mean(mu).data_ttt;
%             
%             figure(mu + (grp-1)*10)
%             subplot(2,1,1)
%             
%             plot(Y')
%             legend(num2str(size(Y,1) - sum( double(isnan( Y(:,1))))))
%             
%             
%             [p, tbl, stats] = friedman(Y',1, 'off');
%             
%             c = multcompare(stats, 'Alpha', alph, 'Display', 'off'); 
%             Test(grp).suj_mean(mu).diff_pourc= [];
% 
%             for j = 1 : size(Y,1)
%                 diff = c( logical( nansum( c(:,1:2) == j, 2)), :);
%                 Test(grp).suj_mean(mu).diff_pourc(j) = length( diff( diff(:,6) < alph))/size(Test(grp).suj_mean(mu).data,1)*100 ;
%             end
% 
%             for i = 1 : length(Y)
%                Test(grp).suj_mean(mu).normality(i) = kstest( Y( :,i));
%                Test(grp).suj_mean(mu).out(:,i) = isoutlier( Y( :,i), 'ThresholdFactor', 6); 
%             end
%             B = [B ; sum( sum( Test(grp).suj_mean(mu).out,2))];
%             
%             %Test(grp).suj_mean(mu).out_pourc_cycl = sum( Test(grp).suj_mean(mu).out,2) / length(Test(grp).suj_mean(mu).out) * 100;%sum( sum( Test(grp).suj_mean(mu).out,2))*100;
%             % Test(grp).suj_mean(mu).out(Test(grp).suj_mean(mu).out_pourc_cycl < 2,:)=0;      % minimum d'outliers pour les compter
%             
%             Test(grp).suj_mean(mu).out_pourc_out = sum(Test(grp).suj_mean(mu).out,2) / sum(sum(Test(grp).suj_mean(mu).out,2))*100;
%             
%             for cycl = 1 : length(Test(grp).suj_mean(mu).out_pourc_out) 
%  
%                 if Test(grp).suj_mean(mu).diff_pourc(cycl) > 50 &&  Test(grp).suj_mean(mu).out_pourc_out(cycl) > 33
%                     Test(grp).suj_mean(mu).outliers(cycl) = 1;
%                 end
% 
%             end   
% 
%             Test(grp).suj_mean(mu).data_ttt( logical( Test(grp).suj_mean(mu).outliers'),:) = NaN;
%         end
%             
%         subplot(2,1,2)
%         plot(Test(grp).suj_mean(mu).data_ttt')
%         legend(num2str(size(Test(grp).suj_mean(mu).data_ttt,1) - sum( double(isnan( Test(grp).suj_mean(mu).data_ttt(:,1))))));
%     end
% end
%% selection de 10 cycles

for grp = 1 : 8
    for suj = 1 : 12 
        for mu = 1 : 9
            
            Data = [];
            RMSE = []; 
            
            Data = All_data(grp).suj(suj).mu(mu).data_norm;
            
            for i = 1 : size(Data,1)
                RMSE(i) = sqrt( mean( ( Data(i,:) - mean(Data) ).^2)); 
            end
            
            [xs, index] = sort(RMSE); 
            
            All_data(grp).suj(suj).mu(mu).data_selec_10 = Data(sort(index(1:10))', :);
        end
    end
end