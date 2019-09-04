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
load([pathname,'\All_data2.mat'])

ordre_mu = [9 8 2 1 3 4 7 5 6];

ordre_mu_name = {'Triceps'; 'Biceps'; 'Anterior Deltoid'; 'Middle Deltoid'; 'Upper Trapezius'; 'Serratus Anterior'; 'Great Pectoral'; 'Extensor Digitorum'; 'Flexor Digitorum' };
ordre_mu_shortname = {'Tri'; 'Bi'; 'AD'; 'MD'; 'UT'; 'SA'; 'GP'; 'EDC'; 'FDS' };
condition_name = {'All-body Struck Staccato'; 'All-body Struck Tenuto'; 'All-body Pressed Staccato'; 'All-body Pressed Tenuto';...
                  'Shoulder Struck Staccato'; 'Shoulder Struck Tenuto'; 'Shoulder Pressed Staccato'; 'Shoulder Pressed Tenuto'};

%%
close all

F=2100/5;
subi = 0;

t = 1/F -630/(2*F) : 1/F : 630/(2*F);

ylim_mu = [0 10; 0 10; 0 4; 0 5; 0 3; 0 2; 0 3; 0 3; 0 3];

%%
for mu = ordre_mu
    subi = subi + 1;
    sujets = [];
    clear('sY','Y','Y1','Y2','Y3','Y4','Yp1','Yp2','Yp3','Yp4','A','B')
    for grp = 1 : 8
        suff = num2str(grp);
%         assignin('base', sprintf('Y%s',suff), All_data(grp).suj_mean(mu).mu_ssAberrant)
%         assignin('base', sprintf('Y%s',suff), All_data(grp).suj_mean(mu).mu_maj)
        assignin('base', sprintf('Y%s',suff), [All_data(grp).suj_mean(mu).mu, [1:12]'])
    end
    
    %% ANOVA3 RM all sujets
    
    Y1(isnan(Y1(:,1)),:) = [];
    Y2(isnan(Y2(:,1)),:) = [];
    Y3(isnan(Y3(:,1)),:) = [];
    Y4(isnan(Y4(:,1)),:) = [];
    Y5(isnan(Y5(:,1)),:) = [];
    Y6(isnan(Y6(:,1)),:) = [];
    Y7(isnan(Y7(:,1)),:) = [];
    Y8(isnan(Y8(:,1)),:) = [];
    
    suj_set = mintersect(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8);

    Y1 = Y1(ismember(Y1(:,end),suj_set),1:end-1);
    Y2 = Y2(ismember(Y2(:,end),suj_set),1:end-1);
    Y3 = Y3(ismember(Y3(:,end),suj_set),1:end-1);
    Y4 = Y4(ismember(Y4(:,end),suj_set),1:end-1);
    Y5 = Y5(ismember(Y5(:,end),suj_set),1:end-1);
    Y6 = Y6(ismember(Y6(:,end),suj_set),1:end-1);
    Y7 = Y7(ismember(Y7(:,end),suj_set),1:end-1);
    Y8 = Y8(ismember(Y8(:,end),suj_set),1:end-1);
    
    
    Y = [Y1; Y2; Y3; Y4; Y5; Y6; Y7; Y8];

    iterations = 1000 ; %1000;
    A = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); ones(size(Y3,1),1); ones(size(Y4,1),1); ...
        zeros(size(Y5,1),1); zeros(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % frappé (0) / pressé (1)
    
    B = [zeros(size(Y1,1),1); ones(size(Y2,1),1); zeros(size(Y3,1),1); ones(size(Y4,1),1); ...
         zeros(size(Y5,1),1); ones(size(Y6,1),1); zeros(size(Y7,1),1); ones(size(Y8,1),1)]; % Staccato (0) / Tenuto (1)
       
    Subj = [suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set]; %Y(:,size(Y,2));
    
    F = spm1d.stats.nonparam.anova2rm(Y, A, B, Subj);
    Fi =F.inference(0.05,'iteration',iterations);
    
    
    %% Plot
    
    Titles = {'A. Struck/Pressed' 'B. Staccato/Tenuto' 'Interaction A/B'};
    
    clu1 = logical([zeros(1,length(suj_set)) zeros(1,length(suj_set)) ones(1,length(suj_set)) ones(1,length(suj_set)) zeros(1,length(suj_set)) zeros(1,length(suj_set)) ones(1,length(suj_set)) ones(1,length(suj_set)); ...
                    zeros(1,length(suj_set)) ones(1,length(suj_set)) zeros(1,length(suj_set)) ones(1,length(suj_set)) zeros(1,length(suj_set)) ones(1,length(suj_set)) zeros(1,length(suj_set)) ones(1,length(suj_set))]');
                
    leg = [{'Pressed' 'Struck'}; {'Tenuto' 'Staccato'}]; 
%     clu_color =[ {'b' 'b' 'c' 'c'}; {'g' 'g' 'r' 'r'}];
    
    clu_color = [0 0 1; 0 1 1; 0 1 0; 0 0.5 0];
    
    for clu = 1 : 2% 7
        figure(1)
        subplot(9,2,subi*2-1+(clu-1))
        hold on

        err_c1 = spm1d.plot.plot_errorcloud(t, mean( Y(clu1(:,clu),:)), std(Y(clu1(:,clu),:)), 'facecolor', clu_color(clu*2-1,:), 'facealpha',0.20) ;
        err_c1.Annotation.LegendInformation.IconDisplayStyle = 'off';
        err_c2 = spm1d.plot.plot_errorcloud(t, mean( Y(~clu1(:,clu),:)), std(Y(~clu1(:,clu),:)), 'facecolor', clu_color(clu*2,:), 'facealpha',0.20);
        err_c2.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
        plot(t, mean( Y(clu1(:,clu),:)), 'Color', clu_color(clu*2-1,:), 'LineWidth', 2.5);
        plot(t, mean( Y(~clu1(:,clu),:)), 'Color', clu_color(clu*2,:), 'LineWidth', 2.5);
        if subi == 1
            legend(leg(clu,:), 'AutoUpdate', 'off')
        end        
%         
%         if clu == 1 
%             ylabel(ordre_mu_shortname{mu},'FontSize',12)
%         end

        xlim([-0.75 0.75])
        ylim(ylim_mu(subi,:))
        if subi~=9
            xticklabels({''})
        end
        
        if subi == 9
            xlabel('Time (s)','FontSize',12)
        end

        delta = 63;
        ind = 1 : size(Fi.SPMs{1, clu}.clusters,2);
        
        i = 1;
        
        while i ~= size(Fi.SPMs{1, clu}.clusters,2)+1 
            
            SP = round(Fi.SPMs{1, clu}.clusters{1, i}.endpoints);
            if SP(1) == 0
                SP(1) = 1;
            end
            
            if SP(1) == SP(2)
                SP(2) = SP(1) + 1;
            end
            
            while i ~= size(Fi.SPMs{1, clu}.clusters,2) && abs(SP(2) - round(Fi.SPMs{1, clu}.clusters{1, i+1}.endpoints(1))) < delta && sum(Fi.SPMs{1, clu}.z(SP(2) : round(Fi.SPMs{1, clu}.clusters{1, i+1}.endpoints(1))) < 3)==0
                    
                SP(2) = round(Fi.SPMs{1, clu}.clusters{1, i+1}.endpoints(2));
                
                if SP(2) == round(Fi.SPMs{1, clu}.clusters{1, size(Fi.SPMs{1, clu}.clusters,2)}.endpoints(2))
                    i = size(Fi.SPMs{1, clu}.clusters,2);
                else
                    i = i+1;
                end
            end
            
            
            if length(SP(1) : SP(2)) > 0.01*630
                hold on
                lim = ylim;
%                 surf(t(SP(1):SP(2)), [lim(1) : (lim(2) - lim(1))/(SP(2)-SP(1)) : lim(2)], repmat(Fi.SPMs{1, clu}.z(SP(1):SP(2)), length(Fi.SPMs{1, clu}.z(SP(1):SP(2))) ,1),'FaceAlpha',0.3, 'EdgeColor','none')
%                 view(2)
                patch([t(SP(1)) t(SP(2)) t(SP(2)) t(SP(1))], [ylim_mu(subi,1), ylim_mu(subi,1), ylim_mu(subi,2), ylim_mu(subi,2)],'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2 )
%                 line(t([SP(1),SP(1)]), ylim,'Color', 'k')
%                 line(t([SP(2),SP(2)]), ylim,'Color', 'k')
%                 line(t([SP(1),SP(2)]), ylim_mu(subi,1)+0.1, 'Color', 'k')
%                 line(t([SP(1),SP(2)]), ylim_mu(subi,2)-0.1, 'Color', 'k')
                
                stats_mu(mu).clu(clu).SP(i,:) = t(SP);
                stats_mu(mu).clu(clu).p = Fi.SPMs{1, clu}.p;  
                stats_mu(mu).clu(clu).d(i) = computeCohen_d( mean(Y(clu1(:,clu),SP(1):SP(2)),2), mean(Y(~clu1(:,clu),SP(1):SP(2)),2) );
                stats_mu(mu).clu(clu).dif(i) = mean(mean(Y(clu1(:,clu),SP(1):SP(2)),2)) - mean(mean(Y(~clu1(:,clu),SP(1):SP(2)),2)) ;
                
                
%                 figure
%                 plot(Y(~clu1(:,clu),SP(1) : SP(2))', 'b')
%                 hold on
%                 plot(Y(clu1(:,clu),SP(1) : SP(2))', 'r')
            end
            
            i = i + 1;
        end
%         subplot(2, 1, 1)
        
%       saveas(figure(clu + 10 * (mu-1)), ['C:\Users\p1218107\Documents\Data_Piano\graphiques\New_reps\Anova3\10 cycles\' All_data(grp).suj_mean(mu).name ' ' leg{clu,1} ' '  leg{clu,2}], 'jpeg')
    end
end

% set(gcf, 'Position', get(0, 'Screensize'))
% saveas(figure(1),['\\10.89.24.15\j\Valentin\Article EMG\figures 2 parameters\ANOVA2'],'jpeg') 