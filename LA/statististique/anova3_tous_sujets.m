%% Initialisation

clear all;
close all;
clc;


%addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))

pathname='C:\Users\p1218107\Documents\Data_Piano\data_newREPS';

%% Ouverture des datas 


cd(pathname)

%load([pathname,'\All_grp.mat'])
load([pathname,'\All_data.mat'])

%%
close all

F=2100/5;

sujet_name={'001','002','003','004','005','006','007','008','009','010','011','012'};


t = 1/F -630/(2*F) : 1/F : 630/(2*F);
for suj = 10 : 12
    for mu = 1 : 9

        sujets = [];
        clear('sY','Y','Y1','Y2','Y3','Y4','Yp1','Yp2','Yp3','Yp4','A','B')
        for grp = 1 : 8
            suff = num2str(grp);

            assignin('base', sprintf('Y%s',suff), All_data(grp).suj(suj).mu(mu).data_selec)

        end

        Y = [Y1; Y2; Y3; Y4; Y5; Y6; Y7; Y8];

        iterations = 1000 ; %1000;
        A = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); ones(size(Y3,1),1); ones(size(Y4,1),1); ...
            zeros(size(Y5,1),1); zeros(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % frappé (0) / pressé (1)

        B = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); zeros(size(Y3,1),1); zeros(size(Y4,1),1); ...
            ones(size(Y5,1),1); ones(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % Bassin (0) / MS (1)

        C = [zeros(size(Y1,1),1); ones(size(Y2,1),1); zeros(size(Y3,1),1); ones(size(Y4,1),1); ...
            zeros(size(Y5,1),1); ones(size(Y6,1),1); zeros(size(Y7,1),1); ones(size(Y8,1),1)]; % Sta (0) / Tenue (1)


        F = spm1d.stats.nonparam.anova3(Y, A, B, C);
        Fi =F.inference(0.01,'iteration',iterations);

    %% preparation data ANOVA3 non répété
    %     
    %     
    %     
    %     Y1 = Y1(randperm(size(Y1,1),7)', 1:end-1);
    %     Y2 = Y2(randperm(size(Y2,1),7)', 1:end-1);
    %     Y3 = Y3(randperm(size(Y3,1),7)', 1:end-1);
    %     Y4 = Y4(randperm(size(Y4,1),7)', 1:end-1);
    %     Y5 = Y5(randperm(size(Y5,1),7)', 1:end-1);
    %     Y6 = Y6(randperm(size(Y6,1),7)', 1:end-1);
    %     Y7 = Y7(randperm(size(Y7,1),7)', 1:end-1);
    %     Y8 = Y8(randperm(size(Y8,1),7)', 1:end-1);
    %    
    %     
    %     Y = [Y1; Y2; Y3; Y4; Y5; Y6; Y7; Y8];
    % 
    %     iterations = 1000 ; %1000;
    %     A = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); ones(size(Y3,1),1); ones(size(Y4,1),1); ...
    %         zeros(size(Y5,1),1); zeros(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % frappé (0) / pressé (1)
    %     
    %     B = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); zeros(size(Y3,1),1); zeros(size(Y4,1),1); ...
    %         ones(size(Y5,1),1); ones(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % Bassin (0) / MS (1)
    %     
    %     C = [zeros(size(Y1,1),1); ones(size(Y2,1),1); zeros(size(Y3,1),1); ones(size(Y4,1),1); ...
    %         zeros(size(Y5,1),1); ones(size(Y6,1),1); zeros(size(Y7,1),1); ones(size(Y8,1),1)]; % Sta (0) / Tenue (1)
    %     
    %     F = spm1d.stats.nonparam.anova3(Y, A, B, C);
    %     Fi =F.inference(0.05,'iteration',iterations);


        %%
        figure(mu)
        set(gcf, 'Position', get(0, 'Screensize'))
        suptitle([All_data(8).suj_mean(mu).name ' suj = ' num2str(suj)])
        subplot(3,3,8:9)
        title('Means All grps')
        hold on
        plot(t, mean(Y1), 'r')
        plot(t, mean(Y2), 'g')
        plot(t, mean(Y3), 'b')
        plot(t, mean(Y4), 'k')
        plot(t, mean(Y5), 'y')
        plot(t, mean(Y6), 'm')
        plot(t, mean(Y7), 'c')
        plot(t, mean(Y8), 'r')


        Titles = {'1. Pressé/Frappé' '2. Bassin / Membre Supérieur' '3. Sta/Ten' 'Interact A/B' 'Interact A/C' 'Interact B/C' 'Interact A/B/C'};

        for clu = 1 : 7
            subplot(3,3,clu)
            plot(t, Fi.SPMs{1, clu}.z,'k')
            hold on 
            line(xlim , [Fi.SPMs{1, clu}.zstar Fi.SPMs{1, clu}.zstar] ,'Color',[1 0 0], 'LineStyle','--')
            title(subplot(3,3,clu), Titles(clu))

            for i = 1 : size(Fi.SPMs{1, clu}.clusters,2)
                SP = round(Fi.SPMs{1, clu}.clusters{1, i}.endpoints);
                if SP(1) == 0
                    SP(1) = 1;
                end
                if SP(2) == 0
                    SP(2) = 1;
                end
                subplot(3,3,8:9)
                hold on
                line(t([SP(1),SP(1)]),ylim,'Color',[1 0 0])
                line(t([SP(2),SP(2)]),ylim,'Color',[1 0 0])
            end
        end

       saveas(figure(mu), ['C:\Users\p1218107\Documents\Data_Piano\graphiques\New_reps\Anova3\par sujet\' sujet_name{suj} '\' All_data(grp).suj_mean(mu).name], 'jpeg')

    end
end

