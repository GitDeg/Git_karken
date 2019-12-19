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
addpath(genpath('\\10.89.24.15\e\Librairies\S2M_Lib\Benjamin_dev\Plots'))
pathname = 'C:\Users\p1218107\Documents\Data_Piano\SAUT';

sujet = {'\001','\002','\003','\004','\005','\006','\007','\008','\009','\010','\011','\012'};

load([pathname '\All_data1s_xy.mat'])

%%
close all
h = figure;
filename = 'SAUT_piano.gif';
%%
data = [];
TF = [];
for i = 1 : 10
    data = [data, All_data(1).CINE(1).mark(5).data_transp(i*3-2:i*3,:)];
    TF = [TF, All_data(1).CINE(1).mark(5).TF(i,:)] ;
%     data = [data, All_data(5).CINE(1).mark(1).data_transp(i*3-2:i*3,:)];
%     TF = [TF, All_data(5).CINE(1).mark(1).TF(i,:)] ;
end
min(data')
max(data')
i = 0;
i_p = 0;
i_i = 0;
 
Z_octave = 0;
Y_octave = [-80 150-80];
l_octave = 170;
Noires = [5.5 -50 Z_octave; 5.5 50 Z_octave; -5.5 50 Z_octave; -5.5 -50 Z_octave];

E4 = [-133.6 -109.3 -109.3 -128.1 -128.1 -133.6; -80 -80 70 70 -30 -30; Z_octave Z_octave Z_octave Z_octave Z_octave Z_octave];
A4 = [-60.71 -36.43 -36.43 -41.93 -41.93 -55.21 -55.21 -60.71; -80 -80 -30 -30 70 70 -30 -30; Z_octave Z_octave Z_octave Z_octave Z_octave Z_octave Z_octave Z_octave];
C5 = [-12.14 12.14 12.14 6.643 6.643 -12.14; -80 -80 -30 -30 70 70; Z_octave Z_octave Z_octave Z_octave Z_octave Z_octave];
E5 = [36.43 60.71 60.71 41.93 41.93 36.43; -80 -80 70 70 -30 -30; Z_octave Z_octave Z_octave Z_octave Z_octave Z_octave];
C6 = [157.9 182.1 182.1 176.6 176.6 157.9; -80 -80 -30 -30 70 70; Z_octave Z_octave Z_octave Z_octave Z_octave Z_octave];

for n = 50:5:length(data)
    % Draw plot for y = x.^n
    x = data(1,n-49:n);
    y = data(2,n-49:n);
    z = data(3,n-49:n);
%     subplot(2,1,1)
    plot3(x,y,z, 'LineWidth', 2) 
    axis equal % this ensures that getframe() returns a consistent size

%% Piano
    hold on 

    for i_piano = 1:3

        X_octave = [-l_octave+l_octave*(i_piano-1) -l_octave+l_octave*i_piano];

        ind = 0;
        for x = X_octave(1) : (X_octave(end)-X_octave(1))/7 : X_octave(end)
            ind = ind+1;
            X_note(ind,:) = [x-(X_octave(end)-X_octave(1))/14 x+(X_octave(end)-X_octave(1))/14];
        end
        X_cadre = [X_note(1,1)  , X_note(1,1)  , X_note(1,1), X_note(end,2) ;X_note(end,2), X_note(end,2), X_note(1,1), X_note(end,2)];
        Y_cadre = [Y_octave, Y_octave(1), Y_octave(1); Y_octave, Y_octave(2), Y_octave(2)];
        Z_cadre = repmat(Z_octave,2,4);

        X_blanches = [X_note(:,1) X_note(:,1)]';
        Y_blanches = [repmat(Y_octave(1), length(X_note),1), repmat(Y_octave(2), length(X_note),1)]';
        Z_blanches = repmat(Z_octave,2, length(X_note));

        X_noires = [Noires(:,1)+X_note(2,1) Noires(:,1)+X_note(3,1) Noires(:,1)+X_note(5,1) Noires(:,1)+X_note(6,1) Noires(:,1)+X_note(7,1)];
        Y_noires = repmat(Noires(:,2)+Y_octave(2)-50,1,5);
        Z_noires = repmat(Noires(:,3),1,5);

        line(X_cadre, Y_cadre, Z_cadre, 'Color', 'k')
        line(X_blanches, Y_blanches, Z_blanches, 'Color', 'k')
        patch(X_noires, Y_noires, Z_noires, 'black')
    end
    
    if sum(TF(n-4:n)) > 0
%         patch(A4(1,:), A4(2,:), A4(3,:), 'b')
        i = i+1;
        if rem(i,2) == 0
            patch(E4(1,:), E4(2,:), E4(3,:), 'y')
            patch(C5(1,:), C5(2,:), C5(3,:), 'y')
        else
            patch(E5(1,:), E5(2,:), E5(3,:), 'y')
            patch(C6(1,:), C6(2,:), C6(3,:), 'y')
        end
    end
    hold off
%%
%     xlim([-90 -10])
%     ylim([-150 70])
%     zlim([0 250])

    xlim([-150 200])
    ylim([-150 70])
    zlim([0 180])
    grid on
    view([-68.505944231614237,12.99878880793003])
    drawnow 
    
%     subplot(2,1,2)
%     axis equal
% %     xlim([All_data(5).CINE(1).mark(1).C(1)-10 All_data(5).CINE(1).mark(1).C(1)+10])
% %     ylim([All_data(5).CINE(1).mark(1).C(2)-30 All_data(5).CINE(1).mark(1).C(2)+30])
%     xlim([All_data(1).CINE(1).mark(5).C1(end,1)-12 All_data(1).CINE(1).mark(5).C2(end,1)+12])
%     ylim([-75 10])
%     if sum(TF(n-4:n)) > 0
%         i = i+1;
%         hold on
% %         plot(All_data(5).CINE(1).mark(1).XY_min(i,1), All_data(5).CINE(1).mark(1).XY_min(i,2), 'rO', 'MarkerSize', 8)
%         if rem(i,2) == 1
%             i_i = i_i+1;
%             plot(All_data(1).CINE(1).mark(5).C2(i_i,1), All_data(1).CINE(1).mark(5).C2(i_i,2), 'rO', 'MarkerSize', 8)
%         else
%             i_p = i_p+1;
%             plot(All_data(1).CINE(1).mark(5).C1(i_i,1), All_data(1).CINE(1).mark(5).C1(i_i,2), 'rO', 'MarkerSize', 8)
%         end
%             
%     end
    
    %% Creation du Gif
    
%  Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
%  Write to the GIF File 
    if n == 50 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end 
end

%%
% close all
% h = figure;
% filename = 'Flexor.gif';
% 
% Fle = [];
% for i = 1 : 5
%     Fle = [Fle, All_data(1).EMG(1).mu(1).data_MVC(i,:)];
% end
% min(Fle');
% max(Fle');
% 
% x = 1 : 5000;
% 
% for n = 50:10:length(Fle)
%     if n < 2000
%         y = Fle(1:n);
%         plot(x(1:n),y, 'LineWidth', 1, 'Color', 'b') 
%         xlim([1 2000])
%         ylim([0 0.40])
%     else
%         plot(x,Fle, 'LineWidth', 1, 'Color', 'b')
%         xlim([n-2000 n])
%         ylim([0 0.40])
%     end
%     grid on
%     drawnow 
% %  Capture the plot as an image 
%     frame = getframe(h); 
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,89); 
% %  Write to the GIF File 
%     if n == 50 
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%     else 
%         imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%     end 
% end

%%
suj = 4
figure(2)
subplot(1,2,1)
t = 0.1 :0.1: 100;
spm1d.plot.plot_meanSD(t,All_data(6).EMG(4).mu(9).data_MVC)
text(10, 0.08, ['CV = ' num2str(All_data(6).EMG(4).mu(9).CV)], 'Fontsize', 12)
xlabel('% cycle')
ylabel('EMG activation (%MVC)')
title('Flexor digitorum superficialis activation')

subplot(1,2,2)
spm1d.plot.plot_meanSD(t,All_data(2).EMG(4).mu(2).data_MVC)
ylim([0 0.035])
text(10, 0.025, ['CV = ' num2str(All_data(2).EMG(4).mu(2).CV)], 'Fontsize', 12)
xlabel('% cycle')
ylabel('EMG activation (%MVC)')
title('Biceps activation')
%%
grpi = 0;

gridx1 = -10:1:10;
gridx2 = -30:2:30;
[x1,x2] = meshgrid(gridx1, gridx2);
xi = [x1(:) x2(:)];

for grp = [1 2 5 6]
    grpi = grpi+1;
    colormap jet
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

        subplot(2,1,1)
%         f1 = ksdensity(C1,xi);
%         surf(x1,x2,reshape(f1,31,21), 'EdgeColor', 'none')
%         contour(x1+79,x2,reshape(f1,31,21))
%         c = colorbar;
%         c.Label.String = 'Probability density function';
%         c.Label.FontSize = 15;
%         caxis([0 0.02])
%         view(2)

        plot(C1(:,1)+79, C1(:,2), 'rO', 'Markersize', 2)
        axis equal
        xlabel('x position (mm)')
        ylabel('y position (mm)')
        
        subplot(2,1,1)
        hold on
%         f2 = ksdensity(C2,xi);
%         surf(x1,x2,reshape(f2,31,21), 'EdgeColor', 'none')
%         contour(x1-79,x2,reshape(f2,31,21))
%         c = colorbar;
%         c.Label.String = 'Probability density function';
%         c.Label.FontSize = 15;
%         caxis([0 0.02])
%         view(2)

        plot(C2(:,1)-79, C2(:,2), 'rO', 'Markersize', 2)
        axis equal
        xlabel('x position (mm)')
        ylabel('y position (mm)')
        xlim([-90 +90])
        ylim([-30 +30])
        text(0,20, ['d= ' num2str(round(mean(ANOVA(10,:,grpi)),1)) '±' num2str(round(std(ANOVA(10,:,grpi)),1))], 'Fontsize', 12,'HorizontalAlignment', 'center')
    else 
        mark = 1;
        C = [];

        for suj = 1 : 12
            C = [C; All_data(grp).CINE(suj).mark(mark).XY_min - All_data(grp).CINE(suj).mark(mark).C]; 
        end
        if grp == 5
            figure(1)
        else
            figure(2)
        end
        
        subplot(2,2,3.5:4)
%         f = ksdensity(C,xi);
%         surf(x1,x2,reshape(f,31,21), 'EdgeColor', 'none')
%         contour(x1,x2,reshape(f,31,21))
%         c = colorbar;
%         c.Label.String = 'Probability density function';
%         c.Label.FontSize = 15;
%         caxis([0 0.02])
%         view(2)
        plot(C(:,1), C(:,2), 'rO', 'Markersize', 2)
        text(0,20, ['d= ' num2str(round(mean(ANOVA(10,:,grpi)),1)) '±' num2str(round(std(ANOVA(10,:,grpi)),1))], 'Fontsize', 12,'HorizontalAlignment', 'center')
        axis equal
        xlim([-10 10])
        ylim([-30 30])
        xlabel('x position (mm)')
        ylabel('y position (mm)')
    end
    
    figure(1)
    suptitle('Endpoint precision without proximal body segment')
    figure(2)
    suptitle('Endpoint precision with proximal body segment')
    
end

%% 

name_parameters = {'Triceps', 'Biceps', 'Anterior deltoid', 'Middle deltoid', 'Upper trapezius', 'Serratus anterior',...
    'Great pectoralis', 'Extensor digitorum communis', 'Flexor digitorum superficialis', 'Precision'};
name_task = {'AH ss bas', 'AH av bas', 'Str ss bas', 'Str av bas'};
ordre_mu = [9 2 1 10];

subi = 0;
for parameter = 3 %ordre_mu % 1:size(ANOVA,1)
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
    
    
    A = [zeros(size(Y1(suj_set),1),1);ones(size(Y2(suj_set),1),1);...
        zeros(size(Y3(suj_set),1),1);ones(size(Y4(suj_set),1),1)];  % 0 = sans bassin   / 1 = avec bassin
    B = [zeros(size(Y1(suj_set),1),1);zeros(size(Y2(suj_set),1),1);...
        ones(size(Y3(suj_set),1),1);ones(size(Y4(suj_set),1),1)];   % 0 = anti-horaire  / 1 = struck
    
    Subj = repmat(suj_set,4,1);
    
    %% Figure
    figure(1)
    if subi < 4
        subplot(3,2,subi*2-1)
        if subi == 3
            boxplot([Y(B==1), Y(B==0)], 'Symbol', 'r', 'Labels', {'Isolated keystrokes', 'Alternated chords'});
            text(1,0.15, ['CV =' num2str(round(mean(Y(B==1)),2)) '±' num2str(round(std(Y(B==1)),2))], 'Fontsize', 12,'HorizontalAlignment', 'center')
            text(2,0.15, ['CV =' num2str(round(mean(Y(B==0)),2)) '±' num2str(round(std(Y(B==0)),2))], 'Fontsize', 12,'HorizontalAlignment', 'center')
        else
            boxplot([Y(B==1), Y(B==0)], 'Symbol', 'r', 'Labels', {'', ''});
            text(1,0.15, ['CV =' num2str(round(mean(Y(B==1)),2)) '±' num2str(round(std(Y(B==1)),2))], 'Fontsize', 12,'HorizontalAlignment', 'center')
            text(2,0.15, ['CV =' num2str(round(mean(Y(B==0)),2)) '±' num2str(round(std(Y(B==0)),2))], 'Fontsize', 12,'HorizontalAlignment', 'center')
        end
        ylim([0 1])
        sigstar([1 2], 0.05) 
        title(name_parameters{parameter})
    else
        subplot(3,2,4)
        boxplot([Y(B==1), Y(B==0)], 'Symbol', 'r', 'Labels', {'Isolated keystrokes', 'Alternated chords'});
        ylim([0 12])
        text(1,1.5, ['d =' num2str(round(mean(Y(B==1)),2)) '±' num2str(round(std(Y(B==1)),2))], 'Fontsize', 12,'HorizontalAlignment', 'center')
        text(2,1.5, ['d =' num2str(round(mean(Y(B==0)),2)) '±' num2str(round(std(Y(B==0)),2))], 'Fontsize', 12,'HorizontalAlignment', 'center')
        sigstar([1 2], 0.05) 
        title(name_parameters{parameter})
    end
    

    %% 
%     figure(3)
%     for i = 1 : 4
%         subplot(2,2,i)
%         boxplot([acti_mu(1).EMG_max(:,:,i); acti_mu(2).EMG_max(:,:,i); acti_mu(3).EMG_max(:,:,i); acti_mu(4).EMG_max(:,:,i); acti_mu(5).EMG_max(:,:,i);...
%              acti_mu(6).EMG_max(:,:,i); acti_mu(7).EMG_max(:,:,i); acti_mu(8).EMG_max(:,:,i); acti_mu(9).EMG_max(:,:,i)]')
%     end
    
end

%% anterior deltoid
boxplot([Y3, Y1, Y4, Y2],  'Symbol', 'r', 'Labels',...
    {'Isolated keystrokes', 'Alternated chords', 'Isolated keystrokes', 'Alternated chords'},...
    'position', [0.8 1.2 1.8 2.2])
ylim([0 1])
sigstar([1.8 2.2], 0.05) 
text(0.8,0.15, ['CV =' num2str(round(nanmean(Y3),2)) '±' num2str(round(nanstd(Y3),2))], 'Fontsize', 12,'HorizontalAlignment', 'center')
text(1.2,0.15, ['CV =' num2str(round(nanmean(Y1),2)) '±' num2str(round(nanstd(Y1),2))], 'Fontsize', 12,'HorizontalAlignment', 'center')
text(1.8,0.15, ['CV =' num2str(round(nanmean(Y4),2)) '±' num2str(round(nanstd(Y4),2))], 'Fontsize', 12,'HorizontalAlignment', 'center')
text(2.2,0.15, ['CV =' num2str(round(nanmean(Y2),2)) '±' num2str(round(nanstd(Y2),2))], 'Fontsize', 12,'HorizontalAlignment', 'center')

title(name_parameters{parameter})