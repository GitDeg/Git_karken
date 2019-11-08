%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                            parameters                             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc 

%% Ubuntu 

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

%% S�lection des Datas � traiter  
%%%%%%%%%%%%% EMG %%%%%%%%%%%%%%

for suj = 1 : 12 
    for grp = 1 : 6
        for mu = 1 : 9    

            Xij = All_data(grp).EMG(suj).mu(mu).data_MVC  ;

        
            %% CV
            All_data(grp).EMG(suj).mu(mu).CV = sqrt(nanmean( nanstd(Xij).*nanstd(Xij) ))/nanmean(nanmean(Xij));
            
            %% EMG max
            
            All_data(grp).EMG(suj).mu(mu).EMG_max = max(Xij');
            
        end
    end
end


%% 
%%%%%%%%% CINE %%%%%%%%%

X = (1:3:30);
Y = (1:3:30) + 1;
Z = (1:3:30) + 2;
x = 1 : 111;
x_rep = repmat(x,10);
grpi = 0;

for grp = [1 2 5 6]
    grpi = grpi+1;
    
    for suj = 1 : 12
        if grp < 5
            mark = 5;
        else 
            mark = 1;
        end
        
        data_X = All_data(grp).CINE(suj).mark(mark).data(X,:);
        data_Y = All_data(grp).CINE(suj).mark(mark).data(Y,:);
        data_Z = All_data(grp).CINE(suj).mark(mark).data(Z,:);
        if (grp == 2 || grp == 1) && suj == 2 
            data_Z = All_data(grp).CINE(suj).mark(mark).data(Z,1:85);
        end

        if grp < 5
            TF = islocalmin(data_Z ,2, 'MinSeparation', 20, 'FlatSelection', 'first', 'MaxNumExtrema', 3);
            TF(:,1) = 1;
        else 
            TF = islocalmin(data_Z ,2, 'FlatSelection', 'first', 'MaxNumExtrema', 1);
        end
        
        All_data(grp).CINE(suj).mark(mark).XY_min = [data_X(TF), data_Y(TF)];

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
           
           figure(1)
           subplot(2,2,grpi)
           hold on
           plot(All_data(grp).CINE(suj).mark(mark).XY_min(:,1), All_data(grp).CINE(suj).mark(mark).XY_min(:,2), '*')
           plot(polyin1)
           plot(x_c1, y_c1, 'r*')
           plot(polyin2)
           plot(x_c2, y_c2, 'r*')
           axis equal
           grid on
           
           All_data(grp).CINE(suj).mark(mark).XY_dist = mean([sqrt((zone1(:,1) - x_c1).^2 + (zone1(:,2) - y_c1).^2) ;...
                                                              sqrt((zone2(:,1) - x_c2).^2 + (zone2(:,2) - y_c2).^2)]);
           
        else
           k = convhull(All_data(grp).CINE(suj).mark(mark).XY_min);
           x = All_data(grp).CINE(suj).mark(mark).XY_min(k,1);
           y = All_data(grp).CINE(suj).mark(mark).XY_min(k,2);
           polyin = polyshape(x,y);
           [x_c , y_c] = centroid(polyin);
           
           figure(1)
           subplot(2,2,grpi)
           hold on
           plot(All_data(grp).CINE(suj).mark(mark).XY_min(:,1), All_data(grp).CINE(suj).mark(mark).XY_min(:,2), '*')
           plot(polyin)
           plot(x_c, y_c, 'r*')
           axis equal
           grid on
           
           All_data(grp).CINE(suj).mark(mark).XY_dist = mean(sqrt((All_data(grp).CINE(suj).mark(mark).XY_min(:,1) - x_c).^2 + (All_data(grp).CINE(suj).mark(mark).XY_min(:,2) - y_c).^2));
        end
    end
end

