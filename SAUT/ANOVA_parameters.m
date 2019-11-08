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

        ANOVA(:,suj,grpi) = [[All_data(grp).EMG(suj).mu.CV], All_data(grp).CINE(suj).mark(mark).XY_dist];

        for mu = 1 : 9
            acti_mu(mu).EMG_max(:,suj,grpi) = mean(All_data(grp).EMG(suj).mu(mu).EMG_max)  ;
        end
    end
end


%%
close all

name_parameters = {'Triceps', 'Biceps', 'DeltAnt', 'DeltMed', 'TrapSup', 'GrandDent',...
    'GrandPec', 'Extenseurs', 'Flechisseurs', 'Precision'};
name_task = {'AH ss bas', 'AH av bas', 'Str ss bas', 'Str av bas'};
ordre_mu = [9 8 2 1 3 4 7 5 6 10];


%% ANOVA 2 sur VR_intra

clear('Y', 'Y1', 'Y2', 'Y3', 'Y4', 'suj_set')
subi = 0;
for parameter = ordre_mu % 1:size(ANOVA,1)
    subi = subi + 1;    
    suj_set = [1:12]';
    
    Y1 = ANOVA(parameter,:,1)';
    Y2 = ANOVA(parameter,:,2)';
    Y3 = ANOVA(parameter,:,3)';
    Y4 = ANOVA(parameter,:,4)';
    
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

    name_parameters{parameter}
    disp_summ(res(parameter).Fi)
    
    %% Figure
    
    figure(1)
    subplot(2,5,subi)
    if subi < 6
        boxplot([Y(A==0), Y(A==1), Y(B==0), Y(B==1)], 'Symbol', 'r', 'PlotStyle', 'compact', 'Labels', {'', '', '', ''});
    else 
        boxplot([Y(A==0), Y(A==1), Y(B==0), Y(B==1)], 'Symbol', 'r', 'PlotStyle', 'compact', 'Labels',...
            name_task);
    end
    
    if subi < 10
        ylim([0 1])
    end
    
    if res(parameter).Fi.SPMs{1,1}.p <= 0.05
    	sigstar([1 2], res(parameter).Fi.SPMs{1,1}.p)
    end
    if res(parameter).Fi.SPMs{1,2}.p <= 0.05
        sigstar([3 4], res(parameter).Fi.SPMs{1,1}.p)
    end
        
    title(name_parameters{parameter})
    
    %%
    figure(2)
    subplot(2,5,subi)
    if subi < 6
        boxplot([Y1, Y2, Y3, Y4], 'Symbol', 'r', 'PlotStyle', 'compact', 'Labels', {'', '', '', ''});
    else 
        boxplot([Y1, Y2, Y3, Y4], 'Symbol', 'r', 'PlotStyle', 'compact', 'Labels',...
            name_task);
    end
    
    if subi < 10
        ylim([0 1])
    end 
    title(name_parameters{parameter})    
    
    %% 
%     figure(3)
%     for i = 1 : 4
%         subplot(2,2,i)
%         boxplot([acti_mu(1).EMG_max(:,:,i); acti_mu(2).EMG_max(:,:,i); acti_mu(3).EMG_max(:,:,i); acti_mu(4).EMG_max(:,:,i); acti_mu(5).EMG_max(:,:,i);...
%              acti_mu(6).EMG_max(:,:,i); acti_mu(7).EMG_max(:,:,i); acti_mu(8).EMG_max(:,:,i); acti_mu(9).EMG_max(:,:,i)]')
%     end
    
end
figure(1)
suptitle('CV et precision')

figure(2)
suptitle('CV et precision')

%%
grpi = 0;

for grp = [1 2 5 6]
    grpi = grpi+1;
    
    for suj = 1 : 12
        if grp < 5
            mark = 5;
        else 
            mark = 1;
        end
        
        if grp < 5
           
           zones = kmeans(All_data(grp).CINE(suj).mark(mark).XY_min,2, 'Replicates', 5);
            
           zone1 = All_data(grp).CINE(suj).mark(mark).XY_min(zones == 1, :);
           k1 = convhull(zone1);
           x1 = zone1(k1,1);
           y1 = zone1(k1,2);
           polyin1 = polyshape(x1 ,y1);
           [x_c1 , y_c1] = centroid(polyin1);
           
           zone2 = All_data(grp).CINE(suj).mark(mark).XY_min(zones == 2, :);
           k2 = convhull(zone2);
           x2 = zone2(k2,1);
           y2 = zone2(k2,2);
           polyin2 = polyshape(x2,y2);
           [x_c2 , y_c2] = centroid(polyin2);
           
           figure(3)
           subplot(2,2,grpi)
           hold on
           plot(All_data(grp).CINE(suj).mark(mark).XY_min(:,1), All_data(grp).CINE(suj).mark(mark).XY_min(:,2), '*')
           plot(polyin1)
           plot(x_c1, y_c1, 'r*')
           plot(polyin2)
           plot(x_c2, y_c2, 'r*')
           axis equal
           grid on
           title(name_task{grpi})
        else
           k = convhull(All_data(grp).CINE(suj).mark(mark).XY_min);
           x = All_data(grp).CINE(suj).mark(mark).XY_min(k,1);
           y = All_data(grp).CINE(suj).mark(mark).XY_min(k,2);
           polyin = polyshape(x,y);
           [x_c , y_c] = centroid(polyin);
           
           figure(3)
           subplot(2,2,grpi)
           hold on
           plot(All_data(grp).CINE(suj).mark(mark).XY_min(:,1), All_data(grp).CINE(suj).mark(mark).XY_min(:,2), '*')
           plot(polyin)
           plot(x_c, y_c, 'r*')
           axis equal
           grid on
           title(name_task{grpi})
        end
    end
end
figure(3)
suptitle('Precision')

%%
t = 1 : 1000;
grpi = 0;
for grp = [1 2 5 6]
    grpi = grpi+1;
    figure(3 + grpi)
    
    if grp < 5
        mark = 5;
    else 
        mark = 1;
    end
        
    for mu = 1 : 9
        subplot(3,3,mu)
        Y = All_data(grp).EMG(1).mean_EMG(mu).data_MVC;
        Y(isnan(Y(:,1)),:)=[];
        spm1d.plot.plot_meanSD(t,Y)
    end
    
end

