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
% 
% for suj = 1 : 12 
%     for grp = 1 : 6
%         for mu = 1 : 9    
% 
%             Xij = All_data(grp).EMG(suj).mu(mu).data_MVC  ;
% 
%         
%             %% CV
%             All_data(grp).EMG(suj).mu(mu).CV = sqrt(nanmean( nanstd(Xij).*nanstd(Xij) ))/nanmean(nanmean(Xij));
%             
%             %% EMG max
%             
%             All_data(grp).EMG(suj).mu(mu).EMG_max = max(Xij');
%             
%         end
%     end
% end


%% 
%%%%%%%%% CINE %%%%%%%%%

X = (1:3:30);
Y = (1:3:30) + 1;
Z = (1:3:30) + 2;
x = 1 : 111;
x_rep = repmat(x,10);
grpi = 0;

% a = rand(1,12);
% b = rand(1,12);
% c = rand(1,12);
cmap = jet(12);
c = 'k';

for grp = [1 2 5 6]
    grpi = grpi+1;
    
    for suj = 1 : 12
        if grp < 5
            mark = 5;
        else 
            mark = 1;
        end

        data_X = All_data(grp).CINE(suj).mark(mark).data_transp(X,:);
        data_Y = All_data(grp).CINE(suj).mark(mark).data_transp(Y,:);
        data_Z = All_data(grp).CINE(suj).mark(mark).data_transp(Z,:);
        
        if (grp == 2 || grp == 1) && suj == 2 
            data_Z = All_data(grp).CINE(suj).mark(mark).data_transp(Z,1:85);
        end

        if grp < 5
            TF = islocalmin(data_Z ,2, 'MinSeparation', 20, 'FlatSelection', 'first', 'MaxNumExtrema', 3);
            TF(:,1) = 1;
        else 
            TF = islocalmin(data_Z ,2, 'FlatSelection', 'first', 'MaxNumExtrema', 1);
        end
        All_data(grp).CINE(suj).mark(mark).TF = TF;
        All_data(grp).CINE(suj).mark(mark).XY_min = [data_X(TF), data_Y(TF)];

        if grp < 5
           
           zones = kmeans(All_data(grp).CINE(suj).mark(mark).XY_min,2, 'Replicates', 5);
            
           zone1 = All_data(grp).CINE(suj).mark(mark).XY_min(zones == 1, :);
           [f1,xi1]=ksdensity(zone1);
           xy1 = xi1(find(max(f1)==f1), :);
           All_data(grp).CINE(suj).mark(mark).C1 = [zone1; xy1];
           
           C1 = nchoosek(1:length(zone1),2);
           for p = 1:length(C1)
               d1(p) = norm(zone1(C1(p,1),:) - zone1(C1(p,2),:));
           end
           
           
           
           zone2 = All_data(grp).CINE(suj).mark(mark).XY_min(zones == 2, :);
           [f2,xi2]=ksdensity(zone2);
           xy2 = xi2(find(max(f2)==f2), :);
           All_data(grp).CINE(suj).mark(mark).C2 = [zone2; xy2]; 
           
           C2 = nchoosek(1:length(zone2),2);
           for p = 1:length(C2)
               d2(p) = norm(zone2(C2(p,1),:) - zone2(C2(p,2),:));
           end
           
%            All_data(grp).CINE(suj).mark(mark).XY_dist = mean([sqrt((zone1(:,1) - xy1(1)).^2 + (zone1(:,2) - xy1(2)).^2) ,...
%                                                               sqrt((zone2(:,1) - xy2(1)).^2 + (zone2(:,2) - xy2(2)).^2)]);
        All_data(grp).CINE(suj).mark(mark).XY_dist = [geomean(d1) ,geomean(d2)];
           
        else
            
           [f,xi]=ksdensity(All_data(grp).CINE(suj).mark(mark).XY_min);
           xy = xi(find(max(f)==f), :);
           All_data(grp).CINE(suj).mark(mark).C = xy;
           
           C = nchoosek(1:length(All_data(grp).CINE(suj).mark(mark).XY_min),2);
           for p = 1:length(C)
               d(p) = norm(All_data(grp).CINE(suj).mark(mark).XY_min(C(p,1),:) - All_data(grp).CINE(suj).mark(mark).XY_min(C(p,2),:));
           end
           
%            All_data(grp).CINE(suj).mark(mark).XY_dist = mean(sqrt((All_data(grp).CINE(suj).mark(mark).XY_min(:,1) - xy(1)).^2 + (All_data(grp).CINE(suj).mark(mark).XY_min(:,2) - xy(2)).^2));
            All_data(grp).CINE(suj).mark(mark).XY_dist = geomean(d);
        end
    end
end

%%
grpi = 0;

gridx1 = -10:1:10;
gridx2 = -30:2:30;
[x1,x2] = meshgrid(gridx1, gridx2);
xi = [x1(:) x2(:)];

for grp = [1 2 5 6]
    grpi = grpi+1;
    
    if grp < 5
        mark = 5;
        C1 = [];
        d1 = [];
        C2 = [];
        d2 = [];

        for suj = 1 : 12
            C1 = [C1; All_data(grp).CINE(suj).mark(mark).C1(1:end-1,:) - All_data(grp).CINE(suj).mark(mark).C1(end,:)];
            d1 = [d1; All_data(grp).CINE(suj).mark(mark).XY_dist(1)];
            C2 = [C2; All_data(grp).CINE(suj).mark(mark).C2(1:end-1,:) - All_data(grp).CINE(suj).mark(mark).C2(end,:)];
            d2 = [d2; All_data(grp).CINE(suj).mark(mark).XY_dist(2)];
        end
        if grp == 1
            figure(1)
        else
            figure(2)
        end
        colormap jet
        subplot(2,2,1)
        f1 = ksdensity(C1,xi);
%         surf(x1,x2,reshape(f1,31,21), 'EdgeColor', 'none')
        contour(x1,x2,reshape(f1,31,21))
        c = colorbar;
        c.Label.String = 'Probability density function';
        c.Label.FontSize = 15;
        caxis([0 0.02])
        view(2)
        axis equal
        xlabel('x position (mm)')
        ylabel('y position (mm)')
        
        subplot(2,2,2)
        f2 = ksdensity(C2,xi);
%         surf(x1,x2,reshape(f2,31,21), 'EdgeColor', 'none')
        contour(x1,x2,reshape(f2,31,21))
        c = colorbar;
        c.Label.String = 'Probability density function';
        c.Label.FontSize = 15;
        caxis([0 0.02])
        view(2)
        axis equal
        xlabel('x position (mm)')
        ylabel('y position (mm)')

    else 
        mark = 1;
        C1 = [];

        for suj = 1 : 12
            C1 = [C1; All_data(5).CINE(1).mark(mark).XY_min - All_data(5).CINE(1).mark(mark).C]; 
        end
        if grp == 5
            figure(1)
        else
            figure(2)
        end
        
        subplot(2,2,3.5:4)
        f = ksdensity(C1,xi);
%         surf(x1,x2,reshape(f,31,21), 'EdgeColor', 'none')
        contour(x1,x2,reshape(f,31,21))
        c = colorbar;
        c.Label.String = 'Probability density function';
        c.Label.FontSize = 15;
        caxis([0 0.02])
        view(2)
        axis equal
        xlabel('x position (mm)')
        ylabel('y position (mm)')
    end
    
    figure(1)
    suptitle('Endpoint precision without proximal body segment')
    figure(2)
    suptitle('Endpoint precision with proximal body segment')
    
end

