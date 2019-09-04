%% Initialisation

clear all;
close all;
clc;

%addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};

GrpNames={'Membre supérieur pressé staccato',...
              'Membre supérieur frappé staccato',...
              'Bassin pressé staccato',...
              'Bassin frappé staccato',...
              'Membre supérieur pressé tenue',...
              'Membre supérieur frappé tenue',...
              'Bassin pressé tenue'...
              'Bassin frappé tenue'};
          
alph = 0.01;


pathname='C:\Users\p1218107\Documents\Data_Piano\data_1.5s\'; 

load([pathname,'Test.mat']);

%%

for suj = 1 : 12
    
    PathName_EMG=[pathname,sujet{suj},'_EMG_GRP.mat'];

    load(PathName_EMG)
    
    for grp = 1 : 8
        for mu = 1 : 10
            
        Test(grp).suj(suj).mu(mu).name = GRP(grp).EMGred(mu).labels;
        Test(grp).suj(suj).mu(mu).data = GRP(grp).EMGred(mu).data; 
        
        end 
    end
end

%%

% 
% for grp = 1:8
%     for suj = 1 : 12
%         for mu = 1 : 10
%             
%             Y = Test(grp).suj(suj).mu(mu).data;
%             
%             [p, tbl, stats] = friedman(Y',1, 'off');
%             c = multcompare(stats, 'Alpha', alpha, 'Display', 'off'); 
%             
%             for j = 1 : size(Y,1)
%                 
%                 diff = c( logical( sum( c(:,1:2) == j, 2)), :);
%                 Test(grp).suj(suj).mu(mu).diff(j) = length( diff( diff(:,6) < alpha)) ;
%             
%             end
%             
%             for i = 1 : length(Y)
%                Test(grp).suj(suj).mu(mu).normality(i) = kstest( Y( :,i));
%                Test(grp).suj(suj).mu(mu).out(:,i) = isoutlier( Y( :,i)); 
%             end
%         end
%     end
% end

%% criteria III 
% 50% different from others and 20% outlying 



% Test(grp).suj(suj).mu(mu).out_pourc = 100*sum(Test(grp).suj(suj).mu(mu).out,2)/sum(sum(Test(grp).suj(suj).mu(mu).out,2));

% for grp = 1 : 8
%     for suj = 1 : 12
%         for mu = 1 : 10 
%             
%             Test(grp).suj(suj).mu(mu).outliers = zeros(1,length(Test(grp).suj(suj).mu(mu).diff_pourc));
%             Test(grp).suj(suj).mu(mu).data_ttt = Test(grp).suj(suj).mu(mu).data;
%             
%             while 1
%                 
%                 Y = Test(grp).suj(suj).mu(mu).data_ttt;
%                 
%                 n1 = size(Y,1) - sum( double(isnan( Y(:,1))));
%             
%                 [p, tbl, stats] = friedman(Y',1, 'off');
%                 c = multcompare(stats, 'Alpha', alpha, 'Display', 'off'); 
% 
%                 for j = 1 : size(Y,1)
% 
%                     diff = c( logical( nansum( c(:,1:2) == j, 2)), :);
%                     Test(grp).suj(suj).mu(mu).diff(j) = length( diff( diff(:,6) < alpha)) ;
% 
%                 end
% 
%                 for i = 1 : length(Y)
%                    Test(grp).suj(suj).mu(mu).normality(i) = kstest( Y( :,i));
%                    Test(grp).suj(suj).mu(mu).out(:,i) = isoutlier( Y( :,i)); 
%                 end
%                 
%                 
%                 for cycl = 1 : length(Test(grp).suj(suj).mu(mu).diff_pourc) 
%                 
%                     if Test(grp).suj(suj).mu(mu).diff_pourc(cycl) > 30 &&  Test(grp).suj(suj).mu(mu).out_pourc(cycl) > 30
%                         Test(grp).suj(suj).mu(mu).outliers(cycl) = 1;
%                     end
% 
%                     if Test(grp).suj(suj).mu(mu).diff_pourc(cycl) > 40 &&  Test(grp).suj(suj).mu(mu).out_pourc(cycl) > 25
%                         Test(grp).suj(suj).mu(mu).outliers(cycl) = 1;
%                     end
% 
%                     if Test(grp).suj(suj).mu(mu).diff_pourc(cycl) > 50 &&  Test(grp).suj(suj).mu(mu).out_pourc(cycl) > 20
%                         Test(grp).suj(suj).mu(mu).outliers(cycl) = 1;
%                     end
% 
%                     if Test(grp).suj(suj).mu(mu).diff_pourc(cycl) > 60 &&  Test(grp).suj(suj).mu(mu).out_pourc(cycl) > 15
%                         Test(grp).suj(suj).mu(mu).outliers(cycl) = 1;
%                     end
% 
%                     if Test(grp).suj(suj).mu(mu).diff_pourc(cycl) > 70 &&  Test(grp).suj(suj).mu(mu).out_pourc(cycl) > 10
%                         Test(grp).suj(suj).mu(mu).outliers(cycl) = 1;
%                     end
% 
%                     if Test(grp).suj(suj).mu(mu).diff_pourc(cycl) > 80 &&  Test(grp).suj(suj).mu(mu).out_pourc(cycl) > 5
%                         Test(grp).suj(suj).mu(mu).outliers(cycl) = 1;
%                     end
% 
%                 end   
%                 
%                 Test(grp).suj(suj).mu(mu).data_ttt( logical( Test(grp).suj(suj).mu(mu).outliers'),:) = NaN;
%                 
%                 n2 = size(Y,1) - sum( double(isnan( Y(:,1))));
%                 
%                 if n2 == n1
%                     break
%                 end
%             end
%             Test(grp).suj(suj).mu(mu).outliers_tot = sum(Test(grp).suj(suj).mu(mu).outliers);
%         end
%     end
% end

%% Normalisation 

for suj = 1 : 12
    for mu = 1 : 10
        
        EMG_mu=[];

        for grp = 1 : 8
           
                EMG_mu=[ EMG_mu; Test(grp).suj(suj).mu(mu).data_ttt];

        end
        
        EMG_mean=nanmean(nanmean(EMG_mu));

        for grp = 1 : 8
            Test(grp).suj(suj).mu(mu).data_norm = Test(grp).suj(suj).mu(mu).data_ttt/EMG_mean;
            Test(grp).suj(suj).mu(mu).mean = nanmean( Test(grp).suj(suj).mu(mu).data_norm); 
        end
    end
end

%% ttt stats sur mean


% Sélection des cas a traiter 
for mu = 6% 1 : 10 
    for grp = 1 : 8 
        
        A = [];
        
        
        for suj = 1 : 12
            A = [A; Test(grp).suj(suj).mu(mu).mean]; 
        end
        
        Test(grp).suj_mean(mu).data = A ; 
        
        Xij = Test(grp).suj_mean(mu).data(:,1:end/3);
        Xi = nanmean(Xij);
        X = nanmean(Xi);

        n_mean = size(Test(grp).suj_mean(mu).data(:,1:788),1); % nb de sujet 

        num = nansum( nansum( (Xij - Xi).^2  / (788*(n_mean-1))));
        den = nansum( nansum( (Xij - X).^2 / (788*n_mean -1)));

        Test(grp).suj_mean(mu).VR_inter = num / den ;
            
        figure(mu)
        subplot(4,2,grp)
        spm1d.plot.plot_meanSD(A)
        hold on 
        plot(A')
        legend(num2str(Test(grp).suj_mean(mu).VR_inter))
        
    end
end
%%
B = [];
close all

for mu = 1 : 10
    for grp = 1 : 8 
  
        Test(grp).suj_mean(mu).outliers = zeros(1,size(Test(grp).suj_mean(mu).data,1));
        Test(grp).suj_mean(mu).data_ttt = Test(grp).suj_mean(mu).data;


            
        if Test(grp).suj_mean(mu).VR_inter > 0.75
        
            Y = Test(grp).suj_mean(mu).data_ttt;
            
            figure(mu + (grp-1)*10)
            subplot(2,1,1)
            
            plot(Y')
            legend(num2str(size(Y,1) - sum( double(isnan( Y(:,1))))))
            
            
            [p, tbl, stats] = friedman(Y',1, 'off');
            
            c = multcompare(stats, 'Alpha', alph, 'Display', 'off'); 
            Test(grp).suj_mean(mu).diff_pourc= [];

            for j = 1 : size(Y,1)
                diff = c( logical( nansum( c(:,1:2) == j, 2)), :);
                Test(grp).suj_mean(mu).diff_pourc(j) = length( diff( diff(:,6) < alph))/size(Test(grp).suj_mean(mu).data,1)*100 ;
            end

            for i = 1 : length(Y)
               Test(grp).suj_mean(mu).normality(i) = kstest( Y( :,i));
               Test(grp).suj_mean(mu).out(:,i) = isoutlier( Y( :,i), 'ThresholdFactor', 6); 
            end
            B = [B ; sum( sum( Test(grp).suj_mean(mu).out,2))];
            
            %Test(grp).suj_mean(mu).out_pourc_cycl = sum( Test(grp).suj_mean(mu).out,2) / length(Test(grp).suj_mean(mu).out) * 100;%sum( sum( Test(grp).suj_mean(mu).out,2))*100;
            % Test(grp).suj_mean(mu).out(Test(grp).suj_mean(mu).out_pourc_cycl < 2,:)=0;      % minimum d'outliers pour les compter
            
            Test(grp).suj_mean(mu).out_pourc_out = sum(Test(grp).suj_mean(mu).out,2) / sum(sum(Test(grp).suj_mean(mu).out,2))*100;
            
            for cycl = 1 : length(Test(grp).suj_mean(mu).out_pourc) 
 
                if Test(grp).suj_mean(mu).diff_pourc(cycl) > 50 &&  Test(grp).suj_mean(mu).out_pourc_out(cycl) > 33
                    Test(grp).suj_mean(mu).outliers(cycl) = 1;
                end

            end   

            Test(grp).suj_mean(mu).data_ttt( logical( Test(grp).suj_mean(mu).outliers'),:) = NaN;
        end
            
        subplot(2,1,2)
        plot(Test(grp).suj_mean(mu).data_ttt')
        legend(num2str(size(Test(grp).suj_mean(mu).data_ttt,1) - sum( double(isnan( Test(grp).suj_mean(mu).data_ttt(:,1))))));
    end
end

%%

B = [];
%close all

mu = 6; grp = 8;

Test(grp).suj_mean(mu).outliers = zeros(1,size(Test(grp).suj_mean(mu).data,1));
Test(grp).suj_mean(mu).data_ttt = Test(grp).suj_mean(mu).data;
    Y = Test(grp).suj_mean(mu).data_ttt(:,1:end);

    figure
    subplot(2,1,1)

    plot(Y')
    legend(num2str(size(Y,1) - sum( double(isnan( Y(:,1))))))


    [p, tbl, stats] = friedman(Y',1, 'off');

    c = multcompare(stats, 'Alpha', alph, 'Display', 'off'); 
    Test(grp).suj_mean(mu).diff_pourc= [];

    for j = 1 : size(Y,1)
        diff = c( logical( nansum( c(:,1:2) == j, 2)), :);
        Test(grp).suj_mean(mu).diff_pourc(j) = length( diff( diff(:,6) < alph))/size(Test(grp).suj_mean(mu).data,1)*100 ;
    end

    for i = 1 : length(Y)
       Test(grp).suj_mean(mu).normality(i) = kstest( Y( :,i));
       Test(grp).suj_mean(mu).out(:,i) = isoutlier( Y( :,i), 'ThresholdFactor', 7); 
    end
    B = [B ; sum( sum( Test(grp).suj_mean(mu).out,2))];

    %Test(grp).suj_mean(mu).out_pourc_cycl = sum( Test(grp).suj_mean(mu).out,2) / length(Test(grp).suj_mean(mu).out) * 100;%sum( sum( Test(grp).suj_mean(mu).out,2))*100;
    % Test(grp).suj_mean(mu).out(Test(grp).suj_mean(mu).out_pourc_cycl < 2,:)=0;      % minimum d'outliers pour les compter

    Test(grp).suj_mean(mu).out_pourc_out = sum(Test(grp).suj_mean(mu).out,2) / sum(sum(Test(grp).suj_mean(mu).out,2))*100;

    for cycl = 1 : 12 

        if Test(grp).suj_mean(mu).diff_pourc(cycl) > 50 &&  Test(grp).suj_mean(mu).out_pourc_out(cycl) >= 40
            Test(grp).suj_mean(mu).outliers(cycl) = 1;
        end

    end   

    Test(grp).suj_mean(mu).data_ttt( logical( Test(grp).suj_mean(mu).outliers'),:) = NaN;


subplot(2,1,2)
plot(Test(grp).suj_mean(mu).data_ttt')
legend(num2str(size(Test(grp).suj_mean(mu).data_ttt,1) - sum( double(isnan( Test(grp).suj_mean(mu).data_ttt(:,1))))));



%%

for grp = 1 : 8
    for mu = 1 : 10
        ind = isnan(Test(grp).suj_mean(mu).data_ttt(:,1));
        SUJ = 1 : 12;
        
        for suj = 1 : 12
            Test(grp).suj(suj).mu(mu).data_ttt2 = Test(grp).suj(suj).mu(mu).data_ttt;
        end
        
        for suj = SUJ(ind)
            Test(grp).suj(suj).mu(mu).data_ttt2 = NaN(1,length(Test(grp).suj(suj).mu(mu).data_ttt)); 
        end
    end
end
%%
for suj = 1 : 12
    for mu = 1 : 10
        
        EMG_mu=[];

        for grp = 1 : 8
                EMG_mu=[ EMG_mu; Test(grp).suj(suj).mu(mu).data_ttt2];
        end
        
        EMG_mean=nanmean(nanmean(EMG_mu));

        for grp = 1 : 8
            Test(grp).suj(suj).mu(mu).data_norm2 = Test(grp).suj(suj).mu(mu).data_ttt2/EMG_mean;
            Test(grp).suj(suj).mu(mu).mean2 = nanmean( Test(grp).suj(suj).mu(mu).data_norm2); 
        end
    end
end

%%
for grp = 1 : 8 
    for mu = 1 : 10
        
        A = [];
        
        for suj = 1 : 12
            if isnan(Test(grp).suj(suj).mu(mu).mean2)
                Test(grp).suj(suj).mu(mu).mean2 = NaN(1,length(Test(grp).suj(suj).mu(mu).mean));
            end
            A = [A; Test(grp).suj(suj).mu(mu).mean2]; 
        end
        
        Post_ttt(grp).suj_mean(mu).data = A ;
        Test(grp).suj_mean(mu).data2 = A ; 
    end
end
%%
for grp = 1 : 8
    Post_ttt(grp).name = GrpNames{grp};
    
    for mu = 1 : 10
        Post_ttt(grp).suj_mean(mu).name = GRP(1).EMG(mu).labels;
    end
end
%%
for grp = 1 : 8
    for suj = 1 : 12
        Post_ttt(grp).suj(suj).name = sujet {suj};
        for mu = 1 : 10
            Post_ttt(grp).suj(suj).mu(mu).name = GRP(1).EMG(mu).labels;
            Post_ttt(grp).suj(suj).mu(mu).data = [Test(grp).suj(suj).mu(mu).data_norm2, ones(size(Test(grp).suj(suj).mu(mu).data_norm2,1),1)*suj]  ;
        end
    end
end

