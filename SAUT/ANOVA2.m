%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                          ANOVA 2 RM                               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc 

%% Ubuntu 
% 
% addpath(genpath('/media/valentin/060C6AFC0C6AE5E1/Users/p1218107/Documents/MATLAB'))
% 
% pathname = '/media/valentin/060C6AFC0C6AE5E1/Users/p1218107/Documents/Data_Piano/SAUT';
% sujet = {'/001','/002','/003','/004','/005','/006','/007','/008','/009','/010','/011','/012'};
% 
% load([pathname '/All_data1s.mat'])

%% Windows
% pathname = '\\10.89.24.15\j\Valentin\SAUT';

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
pathname = 'C:\Users\p1218107\Documents\Data_Piano\SAUT';

sujet = {'\001','\002','\003','\004','\005','\006','\007','\008','\009','\010','\011','\012'};

load([pathname '\All_data1s_xy.mat'])

%% Data preparation

for suj = 1 : 12
    grpi = 0;
    for grp = [1 2 5 6]
            grpi = grpi +1;
            if grp < 5
                mark = 5;
            else
                mark = 1;
            end
            
%             ANOVA_VR_intra(:,suj,grpi) = [[All_data(grp).EMG(suj).mu.VR_intra],...
%                 All_data(grp).CINE(suj).mark(mark).VR_intra_X,...
%                 All_data(grp).CINE(suj).mark(mark).VR_intra_Y,...
%                 All_data(grp).CINE(suj).mark(mark).VR_intra_Z]' ;
%             
            ANOVA_CV(:,suj,grpi) = [All_data(grp).EMG(suj).mu.CV];%,...
%                 All_data(grp).CINE(suj).mark(mark).CV_X,...
%                 All_data(grp).CINE(suj).mark(mark).CV_Y,...
%                 All_data(grp).CINE(suj).mark(mark).CV_Z]' ;
%             
%             for parameter = 1 : 9
%                 ANOVA_CVi(parameter).para(:,suj,grpi) = [All_data(grp).EMG(suj).mu(parameter).CVi];
%             end
% 
%             ANOVA_CVi(10).para(:,suj,grpi) = All_data(grp).CINE(suj).mark(mark).CVi_X;
%             ANOVA_CVi(11).para(:,suj,grpi) = All_data(grp).CINE(suj).mark(mark).CVi_Y;
%             ANOVA_CVi(12).para(:,suj,grpi) = All_data(grp).CINE(suj).mark(mark).CVi_Z;
%             
            ANOVA_XY_dist(suj,grpi) = All_data(grp).CINE(suj).mark(mark).XY_dist  ;
            
            ANOVA(:,suj,grpi) = [[All_data(grp).EMG(suj).mu.CV], All_data(grp).CINE(suj).mark(mark).XY_dist];
    end
end


%%
close all

% F=2100;
% t = 1/F -1000/(2*F) : 1/F : 1000/(2*F);
name_parameters = {'Triceps', 'Biceps', 'DeltAnt', 'DeltMed', 'TrapSup', 'GrandDent',...
    'GrandPec', 'Extenseurs', 'Flechisseurs', 'X axis', 'Y axis', 'Z axis'};
ordre_mu = [9 8 2 1 3 4 7 5 6];


%% ANOVA 2 sur VR_intra

clear('Y', 'Y1', 'Y2', 'Y3', 'Y4', 'suj_set')
for parameter = 1:size(ANOVA_VR_intra,1)
       
    suj_set = [1:12]';
    
    Y1 = ANOVA_VR_intra(parameter,:,1)';
    Y2 = ANOVA_VR_intra(parameter,:,2)';
    Y3 = ANOVA_VR_intra(parameter,:,3)';
    Y4 = ANOVA_VR_intra(parameter,:,4)';
    
    suj_set = suj_set(~isnan(Y1(suj_set)));
    suj_set = suj_set(~isnan(Y2(suj_set)));
    suj_set = suj_set(~isnan(Y3(suj_set)));
    suj_set = suj_set(~isnan(Y4(suj_set)));
    

    Y = [Y1(suj_set); Y2(suj_set); Y3(suj_set); Y4(suj_set)];
    
    iterations = 1000;
    
    A = [zeros(size(Y1(suj_set),1),1);ones(size(Y2(suj_set),1),1);...
        zeros(size(Y3(suj_set),1),1);ones(size(Y4(suj_set),1),1)];  % 0 = sans bassin   / 1 = avec bassin
    B = [zeros(size(Y1(suj_set),1),1);zeros(size(Y2(suj_set),1),1);...
        ones(size(Y3(suj_set),1),1);ones(size(Y4(suj_set),1),1)];   % 0 = anti-horaire  / 1 = struck
    
    Subj = repmat(suj_set,4,1);
    
    F = spm1d.stats.nonparam.anova2rm(Y,A,B,Subj);
    res(parameter).Fi =F.inference(0.05,'iteration',iterations);

    disp_summ(res(parameter).Fi)
    
    %% Figure
    
    figure(1)
    subplot(3,4,parameter)
    boxplot([Y(A==0), Y(A==1), Y(B==0), Y(B==1)], 'Symbol', 'r', 'PlotStyle', 'compact', 'Labels', {'', '', '', ''});
    
    ylim([0 1])
    
    if res(parameter).Fi.SPMs{1,1}.p <= 0.05
    	sigstar([1 2], res(parameter).Fi.SPMs{1,1}.p)
    end
    if res(parameter).Fi.SPMs{1,2}.p <= 0.05
        sigstar([3 4], res(parameter).Fi.SPMs{1,1}.p)
    end
    
    
    title(name_parameters{parameter})
    
end
figure(1)
suptitle('VR intra')

%% ANOVA 2 sur CV
clear('Y', 'Y1', 'Y2', 'Y3', 'Y4', 'suj_set')
for parameter = 1:size(ANOVA_CV,1)
       
    suj_set = [1:12]';

    Y1 = ANOVA_CV(parameter,:,1)';
    Y2 = ANOVA_CV(parameter,:,2)';
    Y3 = ANOVA_CV(parameter,:,3)';
    Y4 = ANOVA_CV(parameter,:,4)';
    
    suj_set = suj_set(~isnan(Y1(suj_set)));
    suj_set = suj_set(~isnan(Y2(suj_set)));
    suj_set = suj_set(~isnan(Y3(suj_set)));
    suj_set = suj_set(~isnan(Y4(suj_set)));
    

    Y = [Y1(suj_set); Y2(suj_set); Y3(suj_set); Y4(suj_set)];
    
    iterations = 1000;
    
    A = [zeros(size(Y1(suj_set),1),1);ones(size(Y2(suj_set),1),1);...
        zeros(size(Y3(suj_set),1),1);ones(size(Y4(suj_set),1),1)];  % 0 = sans bassin   / 1 = avec bassin
    B = [zeros(size(Y1(suj_set),1),1);zeros(size(Y2(suj_set),1),1);...
        ones(size(Y3(suj_set),1),1);ones(size(Y4(suj_set),1),1)];  % 0 = anti-horaire  / 1 = struck
    
    Subj = repmat(suj_set,4,1);
    
    F = spm1d.stats.nonparam.anova2rm(Y,A,B,Subj);
    res(parameter).Fi =F.inference(0.05,'iteration',iterations);

    disp_summ(res(parameter).Fi)
    
    %% Figure
    
    figure(2)
    subplot(3,4,parameter)
    boxplot([Y(A==0), Y(A==1), Y(B==0), Y(B==1)], 'Symbol', 'r', 'PlotStyle', 'compact', 'Labels', {'', '', '', ''});
    
    ylim([0 1])
    
    if res(parameter).Fi.SPMs{1,1}.p <= 0.05
    	sigstar([1 2], res(parameter).Fi.SPMs{1,1}.p)
    end
    if res(parameter).Fi.SPMs{1,2}.p <= 0.05
        sigstar([3 4], res(parameter).Fi.SPMs{1,1}.p)
    end
    
    title(name_parameters{parameter})
end

figure(2)
suptitle('CV')

%% ANOVA 2 sur mean(CVi)
clear('Y', 'Y1', 'Y2', 'Y3', 'Y4', 'suj_set')
for parameter = 1:size(ANOVA_CV,1)
       
    suj_set = [1:12]';

    Y1 = nanmean(ANOVA_CVi(parameter).para(:,:,1))';
    Y2 = nanmean(ANOVA_CVi(parameter).para(:,:,2))';
    Y3 = nanmean(ANOVA_CVi(parameter).para(:,:,3))';
    Y4 = nanmean(ANOVA_CVi(parameter).para(:,:,4))';
    
    suj_set = suj_set(~isnan(Y1(suj_set,1)));
    suj_set = suj_set(~isnan(Y2(suj_set,1)));
    suj_set = suj_set(~isnan(Y3(suj_set,1)));
    suj_set = suj_set(~isnan(Y4(suj_set,1)));
    

    Y = [Y1(suj_set,:); Y2(suj_set,:); Y3(suj_set,:); Y4(suj_set,:)];
    
    iterations = 1000;
    
    A = [zeros(size(Y1(suj_set),1),1);ones(size(Y2(suj_set),1),1);...
        zeros(size(Y3(suj_set),1),1);ones(size(Y4(suj_set),1),1)];  % 0 = sans bassin   / 1 = avec bassin
    B = [zeros(size(Y1(suj_set),1),1);zeros(size(Y2(suj_set),1),1);...
        ones(size(Y3(suj_set),1),1);ones(size(Y4(suj_set),1),1)];  % 0 = anti-horaire  / 1 = struck
    
    Subj = repmat(suj_set,4,1);
    
    F = spm1d.stats.nonparam.anova2rm(Y,A,B,Subj);
    res(parameter).Fi =F.inference(0.05,'iteration',iterations);

    disp_summ(res(parameter).Fi)
    
    %% Figure
    
    figure(3)
    subplot(3,4,parameter)
    boxplot([Y(A==0), Y(A==1), Y(B==0), Y(B==1)], 'Symbol', 'r', 'PlotStyle', 'compact', 'Labels', {'', '', '', ''});
    
    ylim([0 1])
    
    if res(parameter).Fi.SPMs{1,1}.p <= 0.05
    	sigstar([1 2], res(parameter).Fi.SPMs{1,1}.p)
    end
    if res(parameter).Fi.SPMs{1,2}.p <= 0.05
        sigstar([3 4], res(parameter).Fi.SPMs{1,1}.p)
    end
    
    title(name_parameters{parameter})
end
figure(3)
suptitle('mean(CVi)')

% %% ANOVA 2 sur CVi
% clear('Y', 'Y1', 'Y2', 'Y3', 'Y4', 'suj_set')
% for parameter = 1:size(ANOVA_CV,1)
%        
%     suj_set = [1:12]';
% 
%     Y1 = ANOVA_CVi(parameter).para(:,:,1)';
%     Y2 = ANOVA_CVi(parameter).para(:,:,2)';
%     Y3 = ANOVA_CVi(parameter).para(:,:,3)';
%     Y4 = ANOVA_CVi(parameter).para(:,:,4)';
%     
%     suj_set = suj_set(~isnan(Y1(suj_set,1)));
%     suj_set = suj_set(~isnan(Y2(suj_set,1)));
%     suj_set = suj_set(~isnan(Y3(suj_set,1)));
%     suj_set = suj_set(~isnan(Y4(suj_set,1)));
%     
% 
%     Y = [Y1(suj_set,:); Y2(suj_set,:); Y3(suj_set,:); Y4(suj_set,:)];
%     
%     iterations = 1000;
%     
%     A = [zeros(size(Y1(suj_set),1),1);ones(size(Y2(suj_set),1),1);...
%         zeros(size(Y3(suj_set),1),1);ones(size(Y4(suj_set),1),1)];  % 0 = sans bassin   / 1 = avec bassin
%     B = [zeros(size(Y1(suj_set),1),1);zeros(size(Y2(suj_set),1),1);...
%         ones(size(Y3(suj_set),1),1);ones(size(Y4(suj_set),1),1)];  % 0 = anti-horaire  / 1 = struck
%     
%     Subj = repmat(suj_set,4,1);
%     
%     F = spm1d.stats.nonparam.anova2rm(Y,A,B,Subj);
%     res(parameter).Fi =F.inference(0.05,'iteration',iterations);
% 
%     disp_summ(res(parameter).Fi)
%     
%     %% Figure
%     
%     res(parameter).Fi.plot('plot_threshold_label',false, 'plot_p_values',true, 'autoset_ylim',true)
%     suptitle(name_parameters{parameter})
% end
%% ANOVA 2 sur distance XY

clear('Y', 'Y1', 'Y2', 'Y3', 'Y4', 'suj_set')
   
suj_set = [1:12]';

Y1 = ANOVA_XY_dist(:,1);
Y2 = ANOVA_XY_dist(:,2);
Y3 = ANOVA_XY_dist(:,3);
Y4 = ANOVA_XY_dist(:,4);

suj_set = suj_set(~isnan(Y1(suj_set)));
suj_set = suj_set(~isnan(Y2(suj_set)));
suj_set = suj_set(~isnan(Y3(suj_set)));
suj_set = suj_set(~isnan(Y4(suj_set)));


Y = [Y1(suj_set); Y2(suj_set); Y3(suj_set); Y4(suj_set)];

iterations = 1000;

A = [zeros(size(Y1(suj_set),1),1);ones(size(Y2(suj_set),1),1);...
    zeros(size(Y3(suj_set),1),1);ones(size(Y4(suj_set),1),1)];  % 0 = sans bassin   / 1 = avec bassin
B = [zeros(size(Y1(suj_set),1),1);zeros(size(Y2(suj_set),1),1);...
    ones(size(Y3(suj_set),1),1);ones(size(Y4(suj_set),1),1)];  % 0 = anti-horaire  / 1 = struck

Subj = repmat(suj_set,4,1);

F = spm1d.stats.nonparam.anova2rm(Y,A,B,Subj);
Fi =F.inference(0.05,'iteration',iterations);

disp_summ(Fi)

%% Figure

figure(4)
subplot(2,1,1)
boxplot([Y(A==0), Y(A==1), Y(B==0), Y(B==1)], 'Symbol', 'r', 'PlotStyle', 'compact', 'Labels', {'', '', '', ''});
title('distance XY')

subplot(2,1,2)
boxplot([Y1, Y2, Y3, Y4], 'Symbol', 'r', 'PlotStyle', 'compact', 'Labels', {'', '', '', ''});


