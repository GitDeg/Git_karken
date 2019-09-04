%% graphiques variance intra/inter
clear all;
close all;
clc;

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))

load('C:\Users\p1218107\Documents\Data_Piano\data_newREPS\All_data2.mat')

%% Boxplot VR intra_cycle - data

GRP = [];

ordre_mu = [9 8 2 1 3 4 7 5 6];
% ordre_mu_name = {'Flexor Digit'; 'Extensor Digit'; 'Biceps'; 'Triceps'; 'Anterior Deltoid'; 'Middle Deltoid'; 'Great Pectoral'; 'Upper Trapezius'; 'Serratus Anterior'};
%ordre_mu_name = {'Flex Digit'; 'Ext Digit'; 'Bicps'; 'Tricps'; 'Ant Delt'; 'Mid Delt'; 'Grt Pectl'; 'Up Trap'; 'Serr Ant'};
ordre_mu_name = {'FDS'; 'EDC'; 'Bi'; 'Tri'; 'AD'; 'MD'; 'GP'; 'UT'; 'SA'};

condition_name = {'All-body Struck Staccato'; 'All-body Struck Tenuto'; 'All-body Pressed Staccato'; 'All-body Pressed Tenuto';...
                  'Upper-limb Struck Staccato'; 'Upper-limb Struck Tenuto'; 'Upper-limb Pressed Staccato';'Upper-limb Pressed Tenuto'};

for grp = 1 : 8 
    for mu = 1 : 9 
        for suj = 1 : 12
            GRP(grp).boxplot_data(suj,mu) = All_data(grp).suj(suj).mu(ordre_mu(mu)).VR_intra_10  ;
            
            MU(mu).boxplot_data(suj,grp) = All_data(grp).suj(suj).mu(ordre_mu(mu)).VR_intra_10 ;
            
            SUJ(suj).boxplot_data(mu,grp) = All_data(grp).suj(suj).mu(ordre_mu(mu)).VR_intra_10 ;
           
        end
        GRP(grp).boxplot_data_VR_inter(mu) = All_data(grp).suj_mean(ordre_mu(mu)).VR_inter ;
%         MU(mu).boxplot_data_VR_inter_maj(grp) = All_data(grp).suj_mean(mu).VR_inter_maj ;
    end
end
%%
for grp = 1 : 8 
    for mu = 1 : 9
        
        suj_maj = All_data(grp).suj_mean(mu).mu_maj(:,end); 
        
        for suj = 1 : length(suj_maj)
            
            GRP(grp).boxplot_data_maj(suj,mu) = All_data(grp).suj(suj_maj(suj)).mu(mu).VR_intra;
        end
        GRP(grp).boxplot_data_maj(GRP(grp).boxplot_data_maj == 0) = NaN;
    end
end

%% spider plot each suj
close all

for suj = 1 : 12
   
    figure(suj)
    
    radarplot(SUJ(suj).boxplot_data',ordre_mu_name, {'r' 'g' 'b' 'k' 'r' 'g' 'b' 'k'}, {}, {'no','no','no','no',':',':',':',':'},5)
    legend(condition_name, 'Location', 'northeastoutside')
    title(['Sujet ' num2str(suj)])
    set(gcf, 'Position', get(0, 'Screensize'))
%     saveas(figure(suj),['C:\Users\p1218107\Documents\Data_Piano\graphiques\New_reps\Variance inter_intra\',['VR_intra_' num2str(suj)]],'jpg')
    
end

%% spider plot each condition

close all

grp_color = {'r' 'g' 'b' 'k' 'r' 'g' 'b' 'k'};
grp_style = {'no','no','no','no','--','--','--','--'};
grp_style_inter = {'-','-','-','-','--','--','--','--'};


grp_fig = 0;
for grp = 5:8% 8
    mean(GRP(grp).boxplot_data),
    grp_fig = grp_fig +1;
   
    R1 = mean(GRP(grp).boxplot_data) + std(GRP(grp).boxplot_data);
    n1 = size(R1,2);
    R1=[R1 R1(:,1)];
    [Theta1,M]=meshgrid(2*pi/n1*[0:n1]+pi/n1,ones(1,size(R1,1)));
    X1=R1.*sin(Theta1);
    Y1=R1.*cos(Theta1);

    R2 = mean(GRP(grp).boxplot_data) - std(GRP(grp).boxplot_data);
    n2 = size(R2,2);
    R2=[R2 R2(:,1)];
    [Theta2,M]=meshgrid(2*pi/n2*[0:n2]+pi/n2,ones(1,size(R2,1)));
    X2=R2.*sin(Theta2);
    Y2=R2.*cos(Theta2);
    
    subplot(2,2,grp_fig)
    %figure(grp)
    
    hold on;
    F1=fill(X1,Y1,cell2mat(grp_color(grp)),'LineStyle','none');
    F2=fill(X2,Y2,cell2mat({'w'}),'LineStyle','none');
    set(F1,'FaceAlpha',0.2)
    set(F2,'FaceAlpha',1)
    
    radarplot( mean(GRP(grp).boxplot_data),     ordre_mu_name, grp_color(grp), {}, grp_style(grp),5)
%     radarplot( mean(GRP(grp).boxplot_data),     {'' '' '' '' '' '' '' '' ''}, grp_color(grp), {}, grp_style(grp),5)
    
%     inter = GRP(grp).boxplot_data_VR_inter;
%     n_inter=size(inter,2);
%     m_inter=size(inter,1);
%     inter=[inter inter(:,1)];
%     [Theta_inter,M]=meshgrid(2*pi/n_inter*[0:n_inter]+pi/n_inter,ones(1,size(inter,1)));
%     X_inter=inter.*sin(Theta_inter);
%     Y_inter=inter.*cos(Theta_inter);
%     plot(X_inter',Y_inter', grp_color{grp}, 'Marker', 'p', 'LineStyle', grp_style_inter{grp}, 'LineWidth',2 );

%     legend({All_data(:).name}, 'Location', 'northeastoutside')
     title(condition_name(grp))
    
%    saveas(figure(1),['C:\Users\p1218107\Documents\Data_Piano\graphiques\New_reps\Variance inter_intra\',['VR_inter_grp' num2str(grp)]],'jpg')
    
end

%set(gcf, 'Position', get(0, 'Screensize'))
%saveas(figure(1),['\\10.89.24.15\j\Valentin\Article EMG\figures\','spider2'],'jpg')

% subplot(3,3,9)
% 
% radarplot( [1 1 1 1 1 1 1 1 1],ordre_mu_name, {'k'}, {}, {},5)


%%
close all
fig = figure;

grp_color = {'r' 'g' 'b' 'k' 'r' 'g' 'b' 'k'};
grp_style = {'no','no','no','no',':',':',':',':'};
% 
% for grp = 1 : 8
%     boxdata(grp,:) = mean(GRP(grp).boxplot_data);
% end

for grp = 1:8
   
    R1 = mean(GRP(grp).boxplot_data) + std(GRP(grp).boxplot_data);

    R1=[R1 R1(:,1)];
    [Theta,M]=meshgrid(2*pi/n1*[0:n1]+pi/n1,ones(1,size(R1,1)));
    X1=R1.*sin(Theta);
    Y1=R1.*cos(Theta);

    R2 = mean(GRP(grp).boxplot_data) - std(GRP(grp).boxplot_data);

    R2=[R2 R2(:,1)];
    [Theta,M]=meshgrid(2*pi/n1*[0:n1]+pi/n1,ones(1,size(R2,1)));
    X2=R2.*sin(Theta);
    Y2=R2.*cos(Theta);
    
    hold on;
%     F1=fill(X1,Y1,cell2mat(grp_color(grp)),'LineStyle','none');
%     F2=fill(X2,Y2,cell2mat({'w'}),'LineStyle','none');
%     set(F1,'FaceAlpha',0.1)
%     set(F2,'FaceAlpha',1)
    pa = patch([X1 X2], [Y1 Y2], grp_color{grp});
    set(pa,'FaceAlpha',0.15, 'EdgeColor','None' )
    set(get(get(pa,'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');   
end

radarplot(boxdata,ordre_mu_name, grp_color, {}, grp_style,5)
legend(condition_name, 'Location', 'northeastoutside')

set(gcf, 'Position', get(0, 'Screensize'))
%saveas(fig,['C:\Users\p1218107\Documents\Data_Piano\graphiques\New_reps\Variance inter_intra\','VR_intra_all_grp' ],'jpg')

%% figure subplot VR_intra (propre a chaque sujet)

fig1 = figure;

for mu = 1 : 9

    subplot(3,3,mu)
   
    boxplot(MU(mu).boxplot_data)
    ylim([0 1])
    title(All_data(1).suj(1).mu(mu).name)
    
    hold on 
    plot(MU(mu).boxplot_data_VR_inter,'linewidth', 1.5)
    plot(MU(mu).boxplot_data_VR_inter_maj, 'r', 'linewidth', 1.5)
    
    set(gca, 'xtick', [1 2 3 4 5 6 7 8 9])
    set(gca, 'xticklabel', {'BFS';'BFT';'BPS';'BPT';'MFS';'MFT';'MPS';'MPT'})
    

    
%     if mu > 6
%         set(gca, 'xticklabel', {'BaFraSta';'BaFraTen';'BaPreSta';'BaPreTen';'MsFraSta';'MsFraTen';'MsPreSta';'MsPreTen'})
%         set(gca, 'XTickLabelRotation', 30)
%     else 
%         set(gca, 'xticklabel', {'' '' '' '' '' '' '' '' ''})
%     end
end

set(gcf, 'Position', get(0, 'Screensize'))
saveas(fig1,['C:\Users\p1218107\Documents\Data_Piano\graphiques\New_reps\Variance inter_intra\','VR_intra_inter_grp'],'jpg')

%%
fig2 = figure;

for grp = 1 : 8
%     
%     pos = (1:9) - 0.15;
%     pos_maj = (1:9) + 0.15;
    
    subplot(4,2,grp)
    boxplot(GRP(grp).boxplot_data)
    ylim([0 1])
%     boxplot(GRP(grp).boxplot_data,'position', pos ,'PlotStyle','compact', 'Colors', 'b')
%     hold on 
%     boxplot(GRP(grp).boxplot_data_maj,'position', pos_maj ,'PlotStyle','compact', 'Colors', 'r')
    title(All_data(grp).name)
    
    hold on 
    plot([All_data(grp).suj_mean(:).VR_inter],'linewidth', 1.5)
    plot([All_data(grp).suj_mean(:).VR_inter_maj], 'r', 'linewidth', 1.5)
    
    set(gca, 'xtick', [1 2 3 4 5 6 7 8 9 10])
    set(gca, 'xticklabel', {'Tri';'Bi';'DtAnt';'DtMed';'TrpS';'GrdDt';'GrdPc';'Ext';'Fle'})

end
set(gcf, 'Position', get(0, 'Screensize'))
saveas(fig2,['C:\Users\p1218107\Documents\Data_Piano\graphiques\New_reps\Variance inter_intra\','VR_intra_inter_muscle'],'jpg')


%% Figure VR_inter (entre les sujets puis entre les sujets majoritaires)

figure
for grp = 1 : 8
    subplot(4,2,grp)
    plot([All_data(grp).suj_mean(:).VR_inter])
    hold on 
    plot([All_data(grp).suj_mean(:).VR_inter_maj], 'r')
end
%% figure
Colors = [{'r'} {'g'} {'b'} {'k'} {'m'} {'y'} {'c'} {'r'}];
figure 
hold on
for grp = 1 : 8
    pos = (1:9)-0.15+(grp-1)*0.05;
    col = repmat(Colors(grp),1,9);
    boxplot(GRP(grp).boxplot_data,'position', pos ,'PlotStyle','compact', 'Colors', Colors{grp})
end
set(gcf, 'Position', get(0, 'Screensize'))
title('Variation ratio for all tests / each muscle')

set(gca, 'xtick', [1 2 3 4 5 6 7 8 9 10])
set(gca, 'xticklabel', {'Triceps';'Biceps';'DeltAnt';'DeltMed';'TrapSup';'GrandDent';'GrandPec';'Extenseurs';'Flechisseurs'})

%% Complement boxplot
hold on 

for grp = 1 : 8

    plot([sujet(1).grp(grp).mu(1:10).VR_inter], 'Color', Colors{grp}, 'linewidth', 1.5)
    
end 

h = findobj(gca, 'Tag', 'Box');
leg = legend(h([2, 62, 51, 45, 40, 25, 15, 1]),...
             {'Membre supérieur pressé staccato',...
              'Membre supérieur frappé staccato',...
              'Bassin pressé staccato',...
              'Bassin frappé staccato',...
              'Membre supérieur pressé tenue',...
              'Membre supérieur frappé tenue',...
              'Bassin pressé tenue'...
              'Bassin frappé tenue'});
set(leg, 'Position', [0.79 0.2 0 0])

saveas(figure(1),['C:\Users\p1218107\Documents\Data_Piano\graphiques\intervariabilité et synergie\var_inter_intra\','VR_inter_intra_allgrp'],'jpg')

%%

figure
pos = [(1:12)-0.15, (1:12)-0.10, (1:12)+0.10, (1:12)+0.15];
col = [repmat({'r'},1,12) repmat({'g'},1,12) repmat({'b'},1,12) repmat({'k'},1,12)];
boxplot([boxplot_data.grp1',boxplot_data.grp2',boxplot_data.grp3',boxplot_data.grp4'],'position', pos ,'PlotStyle','compact', 'ColorGroup', col)
set(gcf, 'Position', get(0, 'Screensize'))
title('Variation ratio for all tests / each subject')

set(gca, 'xtick', [1 2 3 4 5 6 7 8 9 10 11 12])
set(gca, 'xticklabel', sujet_name)



%% 

for mu = 1 : 10
    for grp = 1 : 8
        for suj = 1 : 12

            GRP(grp).mu(mu).data_gliss(suj,:) = sujet(suj).grp(grp).mu(mu).VR_intra_gliss;
            
        end
        GRP(grp).boxplot_mingliss(:,mu) =  min(GRP(grp).mu(mu).data_gliss')';
    end
end


%%
for suj = 1 : 12
    subplot(3,4, suj); plot(GRP(8).mu(3).data_gliss(suj,:))
    hold on
    ylim([0 1.2])
    line(xlim, [GRP(grp).boxplot_data(suj,3) GRP(grp).boxplot_data(suj,3)], 'Color', 'g', 'Linewidth', 1.5)
    line(xlim, [min(GRP(8).mu(3).data_gliss(suj,:)), min(GRP(8).mu(3).data_gliss(suj,:))], 'Color', 'y', 'Linewidth', 1.5)

end

figure; 
for suj = 1 : 12
    subplot(3, 4, suj);
    plot(sujet(suj).grp(8).mu(3).data(:,1:end-1)')
    hold on 
    plot(mean(sujet(suj).grp(8).mu(3).data(:,1:end-1)), 'r', 'linewidth' , 2)
end
    
%%

close all

for mu = 1:10
    for grp = 1:8
        figure(mu)
        set(gcf, 'Position', get(0, 'Screensize'))
        suptitle(sujet(1).grp(grp).mu(mu).name)
        subplot(4,2,grp)

        hold on 
        plot(GRP(grp).mu(mu).data_gliss', 'k') 
        plot(nanmean(GRP(grp).mu(mu).data_gliss), 'r', 'linewidth', 1.5)
        line(xlim, [nanmean(GRP(grp).boxplot_data(:,mu)) nanmean(GRP(grp).boxplot_data(:,mu))], 'Color', 'g', 'Linewidth', 1.5)
        line(xlim, [nanmean(min(GRP(grp).mu(mu).data_gliss')) nanmean(min(GRP(grp).mu(mu).data_gliss'))], 'Color', 'y', 'Linewidth', 1.5)
        ylim([0 1.2])
    end
    
    saveas(figure(mu), ['C:\Users\p1218107\Documents\Data_Piano\graphiques\intervariabilité et synergie\var_inter_intra\sliding window\' sujet(1).grp(grp).mu(mu).name], 'jpeg')
end
%%
Colors = [{'r'} {'g'} {'b'} {'k'} {'m'} {'y'} {'c'} {'r'}];
figure 
hold on
for grp = 1 : 8
    pos = (1:10)-0.15+(grp-1)*0.05;
    col = repmat(Colors(grp),1,10);
    boxplot(GRP(grp).boxplot_mingliss,'position', pos ,'PlotStyle','compact', 'Colors', Colors{grp})
end
set(gcf, 'Position', get(0, 'Screensize'))
title('Variation ratio for all tests / each muscle')

set(gca, 'xtick', [1 2 3 4 5 6 7 8 9 10])
set(gca, 'xticklabel', {'Fléchisseurs1';'Fléchisseurs2';'Extenseurs';'Biceps';'Triceps';'DeltAnt';'DeltMed';'GrandPec';'TrapSup';'GrandDent'})

%% 
%% nb de synergies optimales : tVAF et VAFm

% data plot tVAF and VAFm
ind = 0;
tVAF = [] ;
VAFm = [] ;
ind_suj = [2:7 9:12]; 
mu = [1 2 4 5 6 9]; 

for suj = ind_suj 
    ind = ind + 1;
    for grp = 1 : 8
        tVAF(grp).data(ind,:)= sujet(suj).grp(grp).tVAF  ; 
    end
end





%% Plot
figure
suptitle('Total variance accounted for (tVAF) and individual muscles variance variations')
set(gcf, 'Position', get(0, 'Screensize'))

for grp = 1 : 8
    subplot(4,2,grp)
    hold on
    spm1d.plot.plot_meanSD(tVAF(grp).data)
    line(xlim,[80,80],'Color','r','LineStyle','--')
    line(xlim,[95,95],'Color','g','LineStyle','-')
    title(GrpNames{grp});
    ylabel('% tVAR')
    xlabel('Number of synergies')
    ylim([50 100])
end

%% Coefficients d'activation tous sujets 
syn_def = 3;

ind = 0;
bar_data= [];

for suj = ind_suj
    ind = ind + 1;
    
    for grp = 1 : 8
        Wgrp = sujet(suj).grp(grp).syn(syn_def).W ./(max(sujet(suj).grp(grp).syn(syn_def).W));
        
        for syn = 1 : syn_def
            bar_data(grp).syn(syn).w(:,ind) = Wgrp(:,syn);
        end
    end
end

%% optimisation des coefficients

n_sujet = size(bar_data(1).syn(1).w,2);

varT = [];

ind_suj = [2:7 9:12];

for grp = 1 : 8

    for extra_loop = 1 : 25

        ni = randperm(n_sujet);

        for i = ni 

            is = ind_suj(i);

            for intra_loop = 1 : 25
                A1 = bar_data(grp).syn(1).w ;
                A2 = bar_data(grp).syn(2).w ;
                A3 = bar_data(grp).syn(3).w ;

%                 H1 = sujet(is).grp(grp).syn(syn_def).H(1,:) ;
%                 H2 = sujet(is).grp(grp).syn(syn_def).H(2,:) ;
%                 H3 = sujet(is).grp(grp).syn(syn_def).H(3,:) ;

                varT(1) = sum( [sum(abs(A1(:,i) - mean(A1(:,ni(ni~=i)),2))) ,... % w1w2w3
                                    sum(abs(A2(:,i) - mean(A2(:,ni(ni~=i)),2))) ,...
                                    sum(abs(A3(:,i) - mean(A3(:,ni(ni~=i)),2)))]);

                varT(2) = sum( [sum(abs(A1(:,i) - mean(A1(:,ni(ni~=i)),2))) ,... % w1w3w2
                                    sum(abs(A3(:,i) - mean(A2(:,ni(ni~=i)),2))) ,...
                                    sum(abs(A2(:,i) - mean(A3(:,ni(ni~=i)),2)))]);

                varT(3) = sum( [sum(abs(A2(:,i) - mean(A1(:,ni(ni~=i)),2))) ,... % w2w1w3
                                    sum(abs(A1(:,i) - mean(A2(:,ni(ni~=i)),2))) ,...
                                    sum(abs(A3(:,i) - mean(A3(:,ni(ni~=i)),2)))]);

                varT(4) = sum( [sum(abs(A2(:,i) - mean(A1(:,ni(ni~=i)),2))) ,... % w2w3w1
                                    sum(abs(A3(:,i) - mean(A2(:,ni(ni~=i)),2))) ,...
                                    sum(abs(A1(:,i) - mean(A3(:,ni(ni~=i)),2)))]);

                varT(5) = sum( [sum(abs(A3(:,i) - mean(A1(:,ni(ni~=i)),2))) ,... % w3w1w2
                                    sum(abs(A1(:,i) - mean(A2(:,ni(ni~=i)),2))) ,...
                                    sum(abs(A2(:,i) - mean(A3(:,ni(ni~=i)),2)))]);

                varT(6) = sum( [sum(abs(A3(:,i) - mean(A1(:,ni(ni~=i)),2))) ,... % w3w2w1
                                    sum(abs(A2(:,i) - mean(A2(:,ni(ni~=i)),2))) ,...
                                    sum(abs(A1(:,i) - mean(A3(:,ni(ni~=i)),2)))]);

                min_var = min(varT);

                if varT(1) == min_var
                    bar_data(grp).syn(1).w(:,i) = A1(:,i);
                    bar_data(grp).syn(2).w(:,i) = A2(:,i);
                    bar_data(grp).syn(3).w(:,i) = A3(:,i);

%                     sujet(is).syn_ord.grp4.H = [H1 ; H2 ; H3];
% 
%                     sujet(is).syn_ord.grp4.W =[A1(:,i) A2(:,i) A3(:,i)] ;


                end

                if varT(2) == min_var
                    bar_data(grp).syn(1).w(:,i) = A1(:,i);
                    bar_data(grp).syn(2).w(:,i) = A3(:,i);
                    bar_data(grp).syn(3).w(:,i) = A2(:,i);

%                     sujet(is).syn_ord.grp4.H = [H1 ; H3 ; H2];
% 
%                     sujet(is).syn_ord.grp4.W =[A1(:,i) A3(:,i) A2(:,i)] ;
                end

                if varT(3) == min_var
                    bar_data(grp).syn(1).w(:,i) = A2(:,i);
                    bar_data(grp).syn(2).w(:,i) = A1(:,i);
                    bar_data(grp).syn(3).w(:,i) = A3(:,i);

%                     sujet(is).syn_ord.grp4.H = [H2 ; H1 ; H3];
% 
%                     sujet(is).syn_ord.grp4.W =[A2(:,i) A1(:,i) A3(:,i)] ;
                end

                if varT(4) == min_var
                    bar_data(grp).syn(1).w(:,i) = A2(:,i);
                    bar_data(grp).syn(2).w(:,i) = A3(:,i);
                    bar_data(grp).syn(3).w(:,i) = A1(:,i);

%                     sujet(is).syn_ord.grp4.H = [H2 ; H3 ; H1];
% 
%                     sujet(is).syn_ord.grp4.W =[A2(:,i) A3(:,i) A1(:,i)] ;
                end

                if varT(5) == min_var
                    bar_data(grp).syn(1).w(:,i) = A3(:,i);
                    bar_data(grp).syn(2).w(:,i) = A1(:,i);
                    bar_data(grp).syn(3).w(:,i) = A2(:,i);

%                     sujet(is).syn_ord.grp4.H = [H3 ; H1 ; H2];
% 
%                     sujet(is).syn_ord.grp4.W =[A3(:,i) A1(:,i) A2(:,i)] ;
                end

                if varT(6) == min_var
                    bar_data(grp).syn(1).w(:,i) = A3(:,i);
                    bar_data(grp).syn(2).w(:,i) = A2(:,i);
                    bar_data(grp).syn(3).w(:,i) = A1(:,i);
 
%                     sujet(is).syn_ord.grp4.H = [H3 ; H2 ; H1];
% 
%                     sujet(is).syn_ord.grp4.W =[A3(:,i) A2(:,i) A1(:,i)] ;
                end
            end

        end
    end    
end

%%

for grp = 1 : 8
    ni = randperm(n_sujet);
    is = ind_suj(i);
    A1 = bar_data(grp).syn(1).w ;
    A2 = bar_data(grp).syn(2).w ;
    A3 = bar_data(grp).syn(3).w ;
    
    if mean(A1(1,:)) > [mean(A2(1,:)) mean(A3(1,:))]
       bar_data(grp).syn(1).w = A1;
       
    elseif mean(A1(6,:)) > [mean(A2(6,:)) mean(A3(6,:))]
        bar_data(grp).syn(3).w = A1;
        
    else 
        bar_data(grp).syn(2).w = A1;
    end
    
    if mean(A2(1,:)) > [mean(A1(1,:)) mean(A3(1,:))]
       bar_data(grp).syn(1).w = A2;
       
    elseif mean(A2(6,:)) > [mean(A1(6,:)) mean(A3(6,:))]
        bar_data(grp).syn(3).w = A2;
        
    else 
        bar_data(grp).syn(2).w = A2;
    end
    
    if mean(A3(1,:)) > [mean(A1(1,:)) mean(A2(1,:))]
       bar_data(grp).syn(1).w = A3;
       
    elseif mean(A3(6,:)) > [mean(A1(6,:)) mean(A2(6,:))]
        bar_data(grp).syn(3).w = A3;
        
    else 
        bar_data(grp).syn(2).w = A3;
    end
end

%% 

for grp = 1 : 8
    for syn = 1 : 3
        bar_data(1).mean_syn(syn).w(:,grp) = mean(bar_data(grp).syn(syn).w')';
    end
end


%% All grp / all participants
% mu = [1 2 4 5 6 9];

for grp = 1 : 8 
    fig = figure(grp);
    suptitle(['Synergy activation coefficients : ' GrpNames{grp}])
    set(gcf, 'Position', get(0, 'Screensize'))
    
    for syn = 1 : syn_def
        subplot(syn_def,1,syn)
        bar(bar_data(grp).syn(syn).w)
        
        set(gca, 'xticklabel',{[]})
        title(['Coefficient synergy ' num2str(syn)])
        ylim([0 1])
    end
    set(gca, 'xticklabel', {'Fléchisseurs1';'Fléchisseurs2';'Biceps';'Triceps';'DeltAnt';'TrapSup'})
    saveas(fig, ['C:\Users\p1218107\Documents\Data_Piano\graphiques\intervariabilité et synergie\synergy\Synergy_grp' num2str(grp)], 'jpeg')
end

%% All grp / mean participant

fig = figure;
suptitle('Synergy activation coefficients : ')
set(gcf, 'Position', get(0, 'Screensize'))

for syn = 1 : syn_def
    subplot(syn_def, 1, syn)
    bar(bar_data(1).mean_syn(syn).w)
    set(gca, 'xticklabel',{[]})
end
set(gca, 'xticklabel', {'Fléchisseurs1';'Fléchisseurs2';'Biceps';'Triceps';'DeltAnt';'TrapSup'})
leg = legend({'Membre supérieur pressé staccato',...
          'Membre supérieur frappé staccato',...
          'Bassin pressé staccato',...
          'Bassin frappé staccato',...
          'Membre supérieur pressé tenue',...
          'Membre supérieur frappé tenue',...
          'Bassin pressé tenue'...
          'Bassin frappé tenue'});
set(leg, 'Position', [0.22 0.54 0 0])

saveas(fig, ['C:\Users\p1218107\Documents\Data_Piano\graphiques\intervariabilité et synergie\synergy\Synergy_all_grp'], 'jpeg')