clear all;
close all;
clc;

%% Ubuntu
% addpath(genpath('/media/valentin/060C6AFC0C6AE5E1/Users/p1218107/Documents/MATLAB'))
% sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};
% pathname='/media/valentin/060C6AFC0C6AE5E1/Users/p1218107/Documents/Data_Piano/LA/MVC/'; 
% cd(pathname)
% load([pathname,'/All_data_LA_1s.mat'])

%% Windows
addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};
pathname='C:\Users\p1218107\Documents\Data_Piano\LA\MVC\'; 
cd(pathname)
load([pathname,'\All_data_LA_1s.mat'])

%% 

ordre_mu = [9 8 2 1 3 4 7 5 6];

ordre_mu_name = {'Triceps'; 'Biceps'; 'Anterior Deltoid'; 'Middle Deltoid'; 'Upper Trapezius'; 'Serratus Anterior'; 'Great Pectoral'; 'Extensor Digitorum'; 'Flexor Digitorum' };
ordre_mu_shortname = {'Tri'; 'Bi'; 'AD'; 'MD'; 'UT'; 'SA'; 'GP'; 'EDC'; 'FDS' };
condition_name = {'All-body Struck Staccato'; 'All-body Struck Tenuto'; 'All-body Pressed Staccato'; 'All-body Pressed Tenuto';...
                  'Shoulder Struck Staccato'; 'Shoulder Struck Tenuto'; 'Shoulder Pressed Staccato'; 'Shoulder Pressed Tenuto'};

%%
close all

F=2100/5;

t = 1/F -630/(2*F) : 1/F : 630/(2*F);

ylim_mu = [0 35; 0 30; 0 5; 0 10; 0 10; 0 7; 0 10; 0 15; 0 15];
subi = 0;

%%
for mu = ordre_mu
    subi = subi + 1;

        sujets = [];
        clear('Y','Y1','Y2','Y3','Y4','Y5','Y6','Y7','Y8')
        
        for grp = 1 : 8
            suff = num2str(grp);
            assignin('base', sprintf('Y%s',suff), [All_data_LA(grp).suj_mean(mu).mEMG*100, [1:12]'])
        end

        if mu == 1 % retrait du sujet 1 pour le triceps car donn�es incoh�rentes 
            Y1(1,1) = NaN;
        end

        %% ANOVA2 RM all sujets

        Y1(isnan(Y1(:,1)),:) = [];
        Y2(isnan(Y2(:,1)),:) = [];
        Y3(isnan(Y3(:,1)),:) = [];
        Y4(isnan(Y4(:,1)),:) = [];
        Y5(isnan(Y5(:,1)),:) = [];
        Y6(isnan(Y6(:,1)),:) = [];
        Y7(isnan(Y7(:,1)),:) = [];
        Y8(isnan(Y8(:,1)),:) = [];

        suj_set_1 = mintersect(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8);

        Y1 = Y1(ismember(Y1(:,end),suj_set_1),1:end-1);
        Y2 = Y2(ismember(Y2(:,end),suj_set_1),1:end-1);
        Y3 = Y3(ismember(Y3(:,end),suj_set_1),1:end-1);
        Y4 = Y4(ismember(Y4(:,end),suj_set_1),1:end-1);
        Y5 = Y5(ismember(Y5(:,end),suj_set_1),1:end-1);
        Y6 = Y6(ismember(Y6(:,end),suj_set_1),1:end-1);
        Y7 = Y7(ismember(Y7(:,end),suj_set_1),1:end-1);
        Y8 = Y8(ismember(Y8(:,end),suj_set_1),1:end-1);


        Y = [Y1; Y2; Y3; Y4; Y5; Y6; Y7; Y8];
       

        iterations = 1500;
        A = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); ones(size(Y3,1),1); ones(size(Y4,1),1); ...
            zeros(size(Y5,1),1); zeros(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % frapp� (0) / press� (1)

        B = [zeros(size(Y1,1),1); ones(size(Y2,1),1); zeros(size(Y3,1),1); ones(size(Y4,1),1); ...
             zeros(size(Y5,1),1); ones(size(Y6,1),1); zeros(size(Y7,1),1); ones(size(Y8,1),1)]; % Staccato (0) / Tenuto (1)

        Subj = [suj_set_1 ;suj_set_1 ;suj_set_1 ;suj_set_1 ;suj_set_1 ;suj_set_1 ;suj_set_1 ;suj_set_1]; %Y(:,size(Y,2));

        F = spm1d.stats.nonparam.anova2rm(Y, A, B, Subj);
        Fi =F.inference(0.05,'iteration',iterations);
        
        p0(1).mu(mu).Fi = Fi;
        p0(1).mu(mu).name = ordre_mu_name{mu};
        
        
        Y_p0d_1 = Y;
        
        %% 
        clear('Y','Y1','Y2','Y3','Y4','Y5','Y6','Y7','Y8')
        
        for grp = 1 : 8
            suff = num2str(grp);
            assignin('base', sprintf('Y%s',suff), [All_data_LA(grp).suj_mean(mu).pEMG*100, [1:12]'])
        end

        if mu == 1 % retrait du sujet 1 pour le triceps car donn�es incoh�rentes 
            Y1(1,1) = NaN;
        end

        
        %% 
        
        Y1(isnan(Y1(:,1)),:) = [];
        Y2(isnan(Y2(:,1)),:) = [];
        Y3(isnan(Y3(:,1)),:) = [];
        Y4(isnan(Y4(:,1)),:) = [];
        Y5(isnan(Y5(:,1)),:) = [];
        Y6(isnan(Y6(:,1)),:) = [];
        Y7(isnan(Y7(:,1)),:) = [];
        Y8(isnan(Y8(:,1)),:) = [];

        suj_set_2 = mintersect(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8);

        Y1 = Y1(ismember(Y1(:,end),suj_set_2),1:end-1);
        Y2 = Y2(ismember(Y2(:,end),suj_set_2),1:end-1);
        Y3 = Y3(ismember(Y3(:,end),suj_set_2),1:end-1);
        Y4 = Y4(ismember(Y4(:,end),suj_set_2),1:end-1);
        Y5 = Y5(ismember(Y5(:,end),suj_set_2),1:end-1);
        Y6 = Y6(ismember(Y6(:,end),suj_set_2),1:end-1);
        Y7 = Y7(ismember(Y7(:,end),suj_set_2),1:end-1);
        Y8 = Y8(ismember(Y8(:,end),suj_set_2),1:end-1);


        Y = [Y1; Y2; Y3; Y4; Y5; Y6; Y7; Y8];
       

        iterations = 1500;
        A = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); ones(size(Y3,1),1); ones(size(Y4,1),1); ...
            zeros(size(Y5,1),1); zeros(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % frapp� (0) / press� (1)

        B = [zeros(size(Y1,1),1); ones(size(Y2,1),1); zeros(size(Y3,1),1); ones(size(Y4,1),1); ...
             zeros(size(Y5,1),1); ones(size(Y6,1),1); zeros(size(Y7,1),1); ones(size(Y8,1),1)]; % Staccato (0) / Tenuto (1)

        Subj = [suj_set_2 ;suj_set_2 ;suj_set_2 ;suj_set_2 ;suj_set_2 ;suj_set_2 ;suj_set_2 ;suj_set_2]; %Y(:,size(Y,2));

        F = spm1d.stats.nonparam.anova2rm(Y, A, B, Subj);
        Fi =F.inference(0.05,'iteration',iterations);
        
        p0(2).mu(mu).Fi = Fi;
        p0(2).mu(mu).name = ordre_mu_name{mu};
        
        Y_p0d_2 = Y;
        
        %% Plot

        Titles = {'Struck', 'Pressed', 'Staccato', 'Tenuto'};

        clu1 = logical([zeros(1,length(suj_set_1)) zeros(1,length(suj_set_1)) ones(1,length(suj_set_1)) ones(1,length(suj_set_1)) zeros(1,length(suj_set_1)) zeros(1,length(suj_set_1)) ones(1,length(suj_set_1)) ones(1,length(suj_set_1)); ...
                        zeros(1,length(suj_set_1)) ones(1,length(suj_set_1)) zeros(1,length(suj_set_1)) ones(1,length(suj_set_1)) zeros(1,length(suj_set_1)) ones(1,length(suj_set_1)) zeros(1,length(suj_set_1)) ones(1,length(suj_set_1))]');
        
        clu2 = logical([zeros(1,length(suj_set_2)) zeros(1,length(suj_set_2)) ones(1,length(suj_set_2)) ones(1,length(suj_set_2)) zeros(1,length(suj_set_2)) zeros(1,length(suj_set_2)) ones(1,length(suj_set_2)) ones(1,length(suj_set_2)); ...
                        zeros(1,length(suj_set_2)) ones(1,length(suj_set_2)) zeros(1,length(suj_set_2)) ones(1,length(suj_set_2)) zeros(1,length(suj_set_2)) ones(1,length(suj_set_2)) zeros(1,length(suj_set_2)) ones(1,length(suj_set_2))]');

        leg = [{'Pressed' 'Struck'}; {'Tenuto' 'Staccato'}]; 
        
        pos = [0.7, 1.3, 2.7, 3.3];
        
        %%
%         for clu = 1 : 2
%             figure(clu)
%             hold on
%             subplot(3,3, subi)
%             
%             hold on
%             yyaxis left
%             boxplot([Y_p0d_1(clu1(:,clu),:), Y_p0d_1(~clu1(:,clu),:), nan(length(Y_p0d_1)/2, 2)], 'Symbol', 'r', 'PlotStyle', 'compact', 'Positions', pos, 'Labels', {'', '', '', ''});
%             
%             a = get(get(gca,'children'),'children');   % Get the handles of all the objects
%             t = get(a,'tag');   % List the names of all the objects 
%             
%             box1 = [a(8); a(12); a(16); a(20); a(24)];   % 
%             box2 = [a(7); a(11); a(15); a(19); a(23)];   % 
%             if clu == 1 % pressed / struck
%                 set(box1, 'Color', 'b');
%                 set(box2, 'Color', 'c');
%             else
%                 set(box1, 'Color', 'g');
%                 set(box2, 'Color', [0 0.5 0]);
%             end
%             
%             yyaxis right
%             boxplot([nan(length(Y_p0d_2)/2, 2), Y_p0d_2(clu2(:,clu),:), Y_p0d_2(~clu2(:,clu),:)], 'Symbol', 'r', 'PlotStyle', 'compact', 'Positions', pos, 'Labels', {'', '', '', ''});
%             
%             a = get(get(gca,'children'),'children');   % Get the handles of all the objects
%             t = get(a,'tag');   % List the names of all the objects 
%              
%             box3 = [a(6); a(10); a(14); a(18); a(22)];
%             box4 = [a(5); a(9); a(13); a(17); a(21)];
%             if clu == 1 % pressed / struck
%                 set(box3, 'Color', 'b');
%                 set(box4, 'Color', 'c');
%             else
%                 set(box3, 'Color', 'g');
%                 set(box4, 'Color', [0 0.5 0]);
%             end
%             
%             if p0(1).mu(mu).Fi.SPMs{1, clu}.p <= 0.05
%                 sigstar(pos(1:2), p0(1).mu(mu).Fi.SPMs{1, clu}.p  )
%             end
%             if p0(2).mu(mu).Fi.SPMs{1, clu}.p <= 0.05
%                 sigstar(pos(3:4), p0(2).mu(mu).Fi.SPMs{1, clu}.p  )
%             end
%             
%             xlabel(ordre_mu_shortname(mu))
%             
%         end
        %%
        l = [0.1 0.29 0.52 0.71];
        b = [0.88 0.77 0.66 0.55 0.44 0.33 0.22 0.11 0.0];
        
        for clu = 1 : 2
            figure(1)
            h = subplot(9,4, subi*4-(4 - (clu*2-1)));
            hold on
            boxplot([Y_p0d_1(clu1(:,clu),:), Y_p0d_1(~clu1(:,clu),:)], 'Symbol', 'r', 'PlotStyle', 'compact', 'Labels', {'', ''}, 'orientation', 'horizontal');
            
            a = get(get(gca,'children'),'children');   % Get the handles of all the objects
            
            box1 = [a(3); a(5); a(7); a(9); a(11)];   % 
            box2 = [a(4); a(6); a(8); a(10); a(12)];   % 
            if clu == 1 % pressed / struck
                set(box1, 'Color', 'b');
                set(box2, 'Color', 'c');
            else
                set(box1, 'Color', 'g');
                set(box2, 'Color', [0 0.5 0]);
            end
            
            xlim([0 round(max([a(12).XData a(11).XData])+ max([a(12).XData a(11).XData])*0.1,0)])
            
            X_lim = max([a(12).XData a(11).XData]) + max([a(12).XData a(11).XData])*0.03;
            
            if p0(1).mu(mu).Fi.SPMs{1, clu}.p <= 0.05
                sigstar_h([1 2], p0(1).mu(mu).Fi.SPMs{1, clu}.p, X_lim)
            end
           
            set(h, 'outerposition', [l(clu*2-1), b(subi), 0.18, 0.1] );
            
            % 
            h = subplot(9,4, subi*4-(4 -clu*2));
            hold on
            boxplot([Y_p0d_2(clu2(:,clu),:), Y_p0d_2(~clu2(:,clu),:)], 'Symbol', 'r', 'PlotStyle', 'compact', 'Labels', {'', ''}, 'orientation', 'horizontal');
            
            a = get(get(gca,'children'),'children');   % Get the handles of all the objects
            box3 = [a(3); a(5); a(7); a(9); a(11)];   % 
            box4 = [a(4); a(6); a(8); a(10); a(12)];
            
            if clu == 1 % pressed / struck
                set(box3, 'Color', 'b');
                set(box4, 'Color', 'c');
            else
                set(box3, 'Color', 'g');
                set(box4, 'Color', [0 0.5 0]);
            end
            
            xlim([0 round(max([a(12).XData a(11).XData])+ max([a(12).XData a(11).XData])*0.1,-1)])
            
            X_lim = max([a(12).XData a(11).XData]) + max([a(12).XData a(11).XData])*0.03;
            
            if p0(2).mu(mu).Fi.SPMs{1, clu}.p <= 0.05
                sigstar([1 2], p0(2).mu(mu).Fi.SPMs{1, clu}.p, X_lim)
            end
            set(h, 'outerposition', [l(clu*2), b(subi), 0.18, 0.1] );
           % xlabel(ordre_mu_shortname(mu))
            
        end
%     suptitle(All_data_LA(grp).suj_mean(mu).name)
end

% set(gcf, 'Position', get(0, 'Screensize'))
% saveas(figure(1),['\\10.89.24.15\j\Valentin\Article EMG\figure MVC\0d\all_mu'],'jpeg') 