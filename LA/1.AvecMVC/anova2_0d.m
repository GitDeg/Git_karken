clear all;
close all;
clc;

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};

pathname='C:\Users\p1218107\Documents\Data_Piano\LA\MVC\'; 

%% Ouverture des datas 

cd(pathname)

%load([pathname,'\All_grp.mat'])
load([pathname,'\All_data_LA_1s.mat'])

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
mu_i = 0;
for mu = ordre_mu
    mu_i = mu_i +1;
    subi = subi + 1;
    for p0d = 1:2
        sujets = [];
        clear('sY','Y','Y1','Y2','Y3','Y4','Yp1','Yp2','Yp3','Yp4','A','B')
        
        if p0d == 1 
            for grp = 1 : 8
                suff = num2str(grp);
                assignin('base', sprintf('Y%s',suff), [All_data_LA(grp).suj_mean(mu).mEMG*100, [1:12]'])
            end

            if mu == 1 % retrait du sujet 1 pour le triceps car données incohérentes 
                Y1(1,1) = NaN;
            end
        else 
            for grp = 1 : 8
                suff = num2str(grp);
                assignin('base', sprintf('Y%s',suff), [All_data_LA(grp).suj_mean(mu).pEMG*100, [1:12]'])
            end

            if mu == 1 % retrait du sujet 1 pour le triceps car données incohérentes 
                Y1(1,1) = NaN;
            end
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

        iterations = 1500;
        A = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); ones(size(Y3,1),1); ones(size(Y4,1),1); ...
            zeros(size(Y5,1),1); zeros(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % frappé (0) / pressé (1)

        B = [zeros(size(Y1,1),1); ones(size(Y2,1),1); zeros(size(Y3,1),1); ones(size(Y4,1),1); ...
             zeros(size(Y5,1),1); ones(size(Y6,1),1); zeros(size(Y7,1),1); ones(size(Y8,1),1)]; % Staccato (0) / Tenuto (1)

        Subj = [suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set]; %Y(:,size(Y,2));

        
        %%
        F = spm1d.stats.nonparam.anova1rm(Y(B==0),A(B==0),Subj(B==0));
        Fi = F.inference(0.05,'iteration',iterations)
        
        %%
        F = spm1d.stats.nonparam.anova2rm(Y, A, B, Subj);
        Fi =F.inference(0.05,'iteration',iterations);
        
        p0(p0d).mu(mu_i).Fi = Fi;
        p0(p0d).mu(mu_i).name = ordre_mu_name{mu};
       
        p0(p0d).mu(mu_i).p(1) = p0(p0d).mu(mu_i).Fi.SPMs{1,1}.p; 
        p0(p0d).mu(mu_i).F(1) = p0(p0d).mu(mu_i).Fi.SPMs{1,1}.z;
        p0(p0d).mu(mu_i).d(1) = computeCohen_d(Y(logical(A==1)), Y(logical(A==0)));
        p0(p0d).mu(mu_i).dif(1) = mean(Y(logical(A==1))) - mean(Y(logical(A==0))) ;
        
        p0(p0d).mu(mu_i).p(2) = p0(p0d).mu(mu_i).Fi.SPMs{1,2}.p;  
        p0(p0d).mu(mu_i).F(2) = p0(p0d).mu(mu_i).Fi.SPMs{1,2}.z;
        p0(p0d).mu(mu_i).d(2) = computeCohen_d(Y(logical(B==1)), Y(logical(B==0)));
        p0(p0d).mu(mu_i).dif(2) = mean(Y(logical(B==1))) - mean(Y(logical(B==0))) ;

        p0(p0d).mu(mu_i).p(3) = p0(p0d).mu(mu_i).Fi.SPMs{1,3}.p;  
        p0(p0d).mu(mu_i).F(3) = p0(p0d).mu(mu_i).Fi.SPMs{1,3}.z;
        %% Plot
% 
%         Titles = {'Struck', 'Pressed', 'Staccato', 'Tenuto'};
% 
%         clu1 = logical([zeros(1,length(suj_set)) zeros(1,length(suj_set)) ones(1,length(suj_set)) ones(1,length(suj_set)) zeros(1,length(suj_set)) zeros(1,length(suj_set)) ones(1,length(suj_set)) ones(1,length(suj_set)); ...
%                         zeros(1,length(suj_set)) ones(1,length(suj_set)) zeros(1,length(suj_set)) ones(1,length(suj_set)) zeros(1,length(suj_set)) ones(1,length(suj_set)) zeros(1,length(suj_set)) ones(1,length(suj_set))]');
% 
%         leg = [{'Pressed' 'Struck'}; {'Tenuto' 'Staccato'}]; 
%         
%         pos = [-0.5, 0.5] + p0d;
% 
%         %%
%         for clu = 1 : 2
%             figure(clu)
%             hold on
%             subplot(3,6, subi*2-(2-p0d))
%             
%             hold on
%             boxplot([Y(clu1(:,clu),:), Y(~clu1(:,clu),:)], 'Symbol', 'r', 'PlotStyle', 'compact');
% 
%             if Fi.SPMs{1, clu}.p <= 0.05
%                 sigstar([1,2], Fi.SPMs{1, clu}.p  )
%             end
%         end
      end
%     suptitle(All_data_LA(grp).suj_mean(mu).name)
end

% set(gcf, 'Position', get(0, 'Screensize'))
% saveas(figure(1),['\\10.89.24.15\j\Valentin\Article EMG\figure MVC\all_mvc_hd_1s2'],'jpeg') 