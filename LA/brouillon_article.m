A = MU(3).boxplot_data(:, [1 3 5 7]); 
B = MU(3).boxplot_data(:, [2 4 6 8]);

%%

mA = mean(A(:))
sA = std(A(:))

mB = mean(B(:))
sB = std(B(:))

[h p] = ttest(A(:), B(:))



%%

% mu = 7

for i_sub = 1 : 9
    for j = 1 : 9
        
%         [h_t(i,j) p_t(i,j)] = ttest(MU(mu).boxplot_data(:, i), MU(mu).boxplot_data(:, j))
        [h_t(i_sub,j) p_t(i_sub,j)] = ttest(mu_data(:, i_sub), mu_data(:, j))
        
    end
end

%%
anova_data = [] ; 

%for mu = 1 : 9
    for grp = 1 : 8
    
    anova_data = [anova_data ; GRP(grp).boxplot_data];
    
    m_grp(grp) = mean(GRP(grp).boxplot_data(:));
    
    
    
    end
    
    %%
%end
B = [zeros(12,1); zeros(12,1); ones(12,1); ones(12,1); ...
    zeros(12,1); zeros(12,1); ones(12,1); ones(12,1)]; % frappé (0) / pressé (1)

A = [zeros(12,1); zeros(12,1); zeros(12,1); zeros(12,1); ...
    ones(12,1); ones(12,1); ones(12,1); ones(12,1)]; % Bassin (0) / MS (1)

C = [zeros(12,1); ones(12,1); zeros(12,1); ones(12,1); ...
    zeros(12,1); ones(12,1); zeros(12,1); ones(12,1)]; % Sta (0) / Tenue (1)

SUBJ = repmat([1:12 ]',8,1);

for mu = 3% : 9
    spm  = spm1d.stats.anova3rm(anova_data(:,mu), A, B, C, SUBJ);
    spmi = spm.inference(0.05);
    disp_summ(spmi)
    
%     for i = 1 : 7
%         MU(mu).table_anova(i).name = spmi.SPMs{1, i}.effect  ;
%         MU(mu).table_anova(i).p_value = spmi.SPMs{1, i}.p  ;
%         if i == 1
%             MU(mu).table_anova(i).mean_sd_0 = [mean(anova_data(A==0,mu)) std(anova_data(A==0,mu))] ;
%             MU(mu).table_anova(i).mean_sd_1 = [mean(anova_data(A==1,mu)) std(anova_data(A==1,mu))] ;
%         end
%         if i == 2
%             MU(mu).table_anova(i).mean_sd_0 = [mean(anova_data(B==0,mu)) std(anova_data(B==0,mu))] ;
%             MU(mu).table_anova(i).mean_sd_1 = [mean(anova_data(B==1,mu)) std(anova_data(B==1,mu))] ;
%         end
%         if i == 3
%             MU(mu).table_anova(i).mean_sd_0 = [mean(anova_data(C==0,mu)) std(anova_data(C==0,mu))] ;
%             MU(mu).table_anova(i).mean_sd_1 = [mean(anova_data(C==1,mu)) std(anova_data(C==1,mu))] ;
%         end
%     end
end

%%

for mu = 1 : 9
    msd(mu,:) = [mean(MU(mu).boxplot_data(:)) std(MU(mu).boxplot_data(:))]
end

%% 
mu_data = [];
for mu = [1 4 2 5 3 6 9 8 7]
    mu_data =  [mu_data MU(mu).boxplot_data(:)];
end
ordre_mu_name = {'Flex Digit'; 'Ext Digit'; 'Bicps'; 'Tricps'; 'Ant Delt'; 'Mid Delt'; 'Grt Pectl'; 'Up Trap'; 'Serr Ant'};
boxplot(mu_data, ordre_mu_name([1 4 2 5 3 6 9 8 7]) )
groups={[1,2],[2,4],[3,4],[4,5],[5,6],[6,8],[7,8],[7,9]};
p_v = [0 0.0000 0.001 0.0106 0.0174 0.0004 0.0021 0.0027];

sigstar(groups, p_v)


%%
mean_value = [0.0719 0.3291 0.4977 0.6838];
grp = [8 11 10 2];
box_data = GRP(7).boxplot_data(:, [1 2 6 8]);  
find(min(abs(box_data - mean_value)) == abs(box_data - mean_value) )
%%

fig = figure 
F=2100/5;
t = 1/F -630/(2*F) : 1/F : 630/(2*F);

subplot(4,1,1)
plot(t, All_data(7).suj(8).mu(9).data_selec_10', 'Linewidth', 1)
title(['Flexor digitorum CV ' 'intra-individual = ' num2str( round( box_data(8,1),3))], 'FontSize',15)


subplot(4,1,2)
plot(t, All_data(7).suj(11).mu(8).data_selec_10', 'Linewidth', 1)
title(['Extensor digitorum CV ' 'intra-individual = ' num2str(round(box_data(11,2),3))], 'FontSize',15)

subplot(4,1,3)
plot(t, All_data(7).suj(10).mu(4).data_selec_10', 'Linewidth', 1)
title(['Middle deltoid CV ' 'intra-individual = ' num2str(round(box_data(10,3),3))], 'FontSize',15)

subplot(4,1,4)
plot(t, All_data(7).suj(2).mu(5).data_selec_10', 'Linewidth', 1)
title(['Upper trapezius CV ' 'intra-individual = ' num2str(round(box_data(2,4),3))], 'FontSize',15)

% set(gcf, 'Position', get(0, 'Screensize'))
% saveas(figure(1),['\\10.89.24.15\j\Valentin\Article EMG\figures\', 'EMG_exemple'],'jpg')

%%

find(min(abs(GRP(7).boxplot_data(:,4)- 0.2256)) == abs(GRP(7).boxplot_data(:,4)- 0.2256))
find(min(abs(GRP(5).boxplot_data(:,4)- 0.6005)) == abs(GRP(5).boxplot_data(:,4)- 0.6005))
 
suj_repr = [10 5];
%%

fig = figure ;
F=2100/5;
t = 1/F -630/(2*F) : 1/F : 630/(2*F);

subplot(2,1,1)
plot(t, All_data(7).suj(10).mu(4).data_selec_10', 'Linewidth', 1.5)
title(['Triceps CV condition 7 ' 'intra-individual = ' num2str( round( GRP(7).boxplot_data(10,4),3))], 'FontSize',20)


subplot(2,1,2)
plot(t, All_data(5).suj(5).mu(4).data_selec_10', 'Linewidth', 1.5)
title(['Triceps CV condition 5 ' 'intra-individual = ' num2str(round( GRP(5).boxplot_data(5,4),3))], 'FontSize',20 )

set(gcf, 'Position', get(0, 'Screensize'))
%saveas(fig,['C:\Users\p1218107\Documents\Data_Piano\graphiques\New_reps\Variance inter_intra\', 'exemple_triceps'],'jpg')

%%

for grp = 1 : 8
    GRP(grp).diff_intra_inter = GRP(grp).boxplot_data_VR_inter - mean(GRP(grp).boxplot_data);
    
end
diff_intra_inter = [] ;
inter = []; 
for grp = 1 : 8
    diff_intra_inter = [diff_intra_inter ; GRP(grp).diff_intra_inter];
    inter = [inter ; GRP(grp).boxplot_data_VR_inter];
end

for i_sub = 1 : 9
    for j = 1 : 9
        
%         [h_t(i,j) p_t(i,j)] = ttest(MU(mu).boxplot_data(:, i), MU(mu).boxplot_data(:, j))
        [h_t(i_sub,j) p_t(i_sub,j)] = ttest(diff_intra_inter(:, i_sub), diff_intra_inter(:, j))
        
    end
end

%%


figure 

sujets = [];
clear('sY','Y','Y1','Y2','Y3','Y4','Yp1','Yp2','Yp3','Yp4','A','B')
for grp = 1 : 8
    suff = num2str(grp);

%         assignin('base', sprintf('Y%s',suff), All_data(grp).suj_mean(mu).mu_ssAberrant)
%         assignin('base', sprintf('Y%s',suff), All_data(grp).suj_mean(mu).mu_maj)
    assignin('base', sprintf('Y%s',suff), [All_data(grp).suj_mean(mu).mu, [1:12]'])
end

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


% Y = [Y1; Y2; Y3; Y4; Y5; Y6; Y7; Y8];

Y = [Y2;Y4;Y6;Y8];

iterations = 1000 ; %1000;

A = [zeros(12,1); ones(12,1); zeros(12,1); ones(12,1)];
Subj = [suj_set ;suj_set ;suj_set ;suj_set];

spm       = spm1d.stats.anova1rm(Y, A, Subj);  %within-subjects model
spmi      = spm.inference(0.05);
disp(spmi)



%(2) Plot:
subplot(2,1,2)
spmi.plot();
spmi.plot_threshold_label();
spmi.plot_p_values();
hold on
plot(spm_bs.z, 'r')  %between-subjects model

%%
close all
F=2100/5;
t = 1/F -630/(2*F) : 1/F : 630/(2*F);
suj_set = [1:12]';
mu_i = 0;

%%

for mu = [9 8 2 1 3 4 7 5 6] %[2 3 4 7 5] % [9 8 2 1 5] % 5 ;
    clu = 1
    mu_i = mu_i+1;

    for grp = 1 : 8
        suff = num2str(grp);

    %         assignin('base', sprintf('Y%s',suff), All_data(grp).suj_mean(mu).mu_ssAberrant)
    %         assignin('base', sprintf('Y%s',suff), All_data(grp).suj_mean(mu).mu_maj)
        assignin('base', sprintf('Y%s',suff), [All_data(grp).suj_mean(mu).mu])
    end

    Y = [Y1; Y2; Y3; Y4; Y5; Y6; Y7; Y8];

    A = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); ones(size(Y3,1),1); ones(size(Y4,1),1); ...
            zeros(size(Y5,1),1); zeros(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % frappé (0) / pressé (1)

    B = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); zeros(size(Y3,1),1); zeros(size(Y4,1),1); ...
        ones(size(Y5,1),1); ones(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % Bassin (0) / MS (1)

    C = [zeros(size(Y1,1),1); ones(size(Y2,1),1); zeros(size(Y3,1),1); ones(size(Y4,1),1); ...
        zeros(size(Y5,1),1); ones(size(Y6,1),1); zeros(size(Y7,1),1); ones(size(Y8,1),1)]; % Sta (0) / Tenue (1)

    Subj = [suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set]; %Y(:,size(Y,2));

    F = spm1d.stats.nonparam.anova3rm(Y, A, B, C, Subj);
    Fi =F.inference(0.05,'iteration',1000);

    fig = figure(1);
    subplot(9,3,mu_i*3)
    hold on

    clu1 = logical([zeros(1,length(suj_set)) zeros(1,length(suj_set)) ones(1,length(suj_set)) ones(1,length(suj_set)) zeros(1,length(suj_set)) zeros(1,length(suj_set)) ones(1,length(suj_set)) ones(1,length(suj_set)); ...
                    zeros(1,length(suj_set)) zeros(1,length(suj_set)) zeros(1,length(suj_set)) zeros(1,length(suj_set)) ones(1,length(suj_set)) ones(1,length(suj_set)) ones(1,length(suj_set)) ones(1,length(suj_set)); ...
                    zeros(1,length(suj_set)) ones(1,length(suj_set)) zeros(1,length(suj_set)) ones(1,length(suj_set)) zeros(1,length(suj_set)) ones(1,length(suj_set)) zeros(1,length(suj_set)) ones(1,length(suj_set))]');

    leg = [{'Pressed' 'Struck'}; {'Shoulder' 'All-body'} ;{'Tenuto' 'Staccato'}]; 
    clu_color =[ {'b' 'k' 'r' 'k'}; {':r' 'k' 'r' 'k' }; {'g' 'k' 'r' 'k'}];


     

    
    err_c1 = spm1d.plot.plot_errorcloud(t, mean( Y(clu1(:,clu),:)), std(Y(clu1(:,clu),:)), 'facecolor', clu_color{clu,2}, 'facealpha',0.20) ;
    err_c1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    err_c2 = spm1d.plot.plot_errorcloud(t, mean( Y(~clu1(:,clu),:)), std(Y(~clu1(:,clu),:)), 'facecolor', clu_color{clu,4}, 'facealpha',0.20);
    err_c2.Annotation.LegendInformation.IconDisplayStyle = 'off';


    plot(t, mean( Y(clu1(:,clu),:)), clu_color{clu,1},'LineWidth', 3);
    plot(t, mean( Y(~clu1(:,clu),:)), clu_color{clu,3},'LineWidth', 3);
    %legend(leg(clu,:), 'AutoUpdate', 'off')
    yyaxis left
    %ylabel(['EMG activation'], 'FontSize',11)
    %ylabel([All_data(1).suj_mean(mu).name ' activation'], 'FontSize',15)

    yyaxis right
    plot(t, Fi.SPMs{1, clu}.z)%, 'color', [0 0 0 0.5])
    hold on 
    line(xlim , [Fi.SPMs{1, clu}.zstar Fi.SPMs{1, clu}.zstar] ,'Color',[1 0 0 0.5], 'LineStyle','--')
    %ylabel('SPM \{F\}', 'FontSize',15)


    delta = 78;
    ind = 1 : size(Fi.SPMs{1, clu}.clusters,2);

    i_sub = 1;

    while i_sub ~= size(Fi.SPMs{1, clu}.clusters,2)+1 

        SP = round(Fi.SPMs{1, clu}.clusters{1, i_sub}.endpoints);
        if SP(1) == 0
            SP(1) = 1;
        end

        if SP(1) == SP(2)
            SP(2) = SP(1) + 1;
        end

        while i_sub ~= size(Fi.SPMs{1, clu}.clusters,2) && abs(SP(2) - round(Fi.SPMs{1, clu}.clusters{1, i_sub+1}.endpoints(1))) < delta
            SP(2) = round(Fi.SPMs{1, clu}.clusters{1, i_sub+1}.endpoints(2));

            if SP(2) == round(Fi.SPMs{1, clu}.clusters{1, size(Fi.SPMs{1, clu}.clusters,2)}.endpoints(2))
                i_sub = size(Fi.SPMs{1, clu}.clusters,2);
            else
                i_sub = i_sub+1;
            end
        end

        hold on
        lim = ylim;
        line(t([SP(1),SP(1)]),ylim,'Color',[0 0 0], 'Linewidth', 1.5)
        line(t([SP(2),SP(2)]),ylim,'Color',[0 0 0], 'Linewidth', 1.5)
        i_sub = i_sub + 1;
    end

end
%xlabel('')
%xlabel('Time (s)', 'Fontsize', 15)

%% Interaction 

% FDS interaction touch / articulation
% Tri interaction touch / body_implication

mu = 1 ;

sujets = [];
clear('sY','Y','Y1','Y2','Y3','Y4','Yp1','Yp2','Yp3','Yp4','A','B')
for grp = 1 : 8
    suff = num2str(grp);

%         assignin('base', sprintf('Y%s',suff), All_data(grp).suj_mean(mu).mu_ssAberrant)
%         assignin('base', sprintf('Y%s',suff), All_data(grp).suj_mean(mu).mu_maj)
    assignin('base', sprintf('Y%s',suff), [All_data(grp).suj_mean(mu).mu])
end

Y = [Y1; Y2; Y3; Y4; Y5; Y6; Y7; Y8];

A = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); ones(size(Y3,1),1); ones(size(Y4,1),1); ...
    zeros(size(Y5,1),1); zeros(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % frappé (0) / pressé (1)

B = [zeros(size(Y1,1),1); zeros(size(Y2,1),1); zeros(size(Y3,1),1); zeros(size(Y4,1),1); ...
    ones(size(Y5,1),1); ones(size(Y6,1),1); ones(size(Y7,1),1); ones(size(Y8,1),1)]; % Bassin (0) / MS (1)

C = [zeros(size(Y1,1),1); ones(size(Y2,1),1); zeros(size(Y3,1),1); ones(size(Y4,1),1); ...
    zeros(size(Y5,1),1); ones(size(Y6,1),1); zeros(size(Y7,1),1); ones(size(Y8,1),1)]; % Sta (0) / Tenue (1)

suj_set = [1:12]';

Subj = [suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set]; %Y(:,size(Y,2));

F = spm1d.stats.nonparam.anova1rm(Y(A==0,:), B(A==0), Subj(A==0));
Fi =F.inference(0.05,'iteration',1000);

figure(2)
subplot(2,1,1)
Fi.plot();

%% interaction CV

% MD biceps touch x articulation

% SA Body x articulation
%end
A = [zeros(12,1); zeros(12,1); ones(12,1); ones(12,1); ...
    zeros(12,1); zeros(12,1); ones(12,1); ones(12,1)]; % frappé (0) / pressé (1)

B = [zeros(12,1); zeros(12,1); zeros(12,1); zeros(12,1); ...
    ones(12,1); ones(12,1); ones(12,1); ones(12,1)]; % Bassin (0) / MS (1)

C = [zeros(12,1); ones(12,1); zeros(12,1); ones(12,1); ...
    zeros(12,1); ones(12,1); zeros(12,1); ones(12,1)]; % Sta (0) / Tenue (1)

SUBJ = repmat([1:12 ]',8,1);

mu = 9;
spm  = spm1d.stats.anova1rm(anova_data(A==1,mu), B(A==1), SUBJ(A==1));
spmi = spm.inference(0.05);
disp(spmi)
    
%% visualisation variabilité 
% Touch

mu = 5;
ordre_mu = [9 8 2 1 3 4 7 5 6];
F=2100/5;
t = 1/F -630/(2*F) : 1/F : 630/(2*F);

m_p = mean( mean( MU(mu).boxplot_data(:,[3 4 7 8])));
m_s = mean( mean( MU(mu).boxplot_data(:,[1 2 5 6])));

[a b] = min(abs(MU(mu).boxplot_data(:, [3 4 7 8]) - m_p) + abs(MU(mu).boxplot_data(:, [1 2 5 6]) - m_s));
[c d] = min(a);
suj = b(d);
grp = [3 4 7 8 ; 1 2 5 6];

figure 
subplot(2,1,1)
hold on
title('pressed')
plot(t,All_data(grp(1,d)).suj(suj).mu(ordre_mu(mu)).data_selec_10')

subplot(2,1,2)
hold on
title('struck')
plot(t,All_data(grp(2,d)).suj(suj).mu(ordre_mu(mu)).data_selec_10')

%% visualisation variabilité 
% Articulation

mu = 5;
ordre_mu = [9 8 2 1 3 4 7 5 6];
F=2100/5;
t = 1/F -630/(2*F) : 1/F : 630/(2*F);

m_s = mean( mean( MU(mu).boxplot_data(:,[1 3 5 7])));
m_t = mean( mean( MU(mu).boxplot_data(:,[2 4 6 8])));

[a b] = min(abs(MU(mu).boxplot_data(:, [1 3 5 7]) - m_s) + abs(MU(mu).boxplot_data(:, [2 4 6 8]) - m_t));
[c d] = min(a);
suj = b(d);
grp = [1 3 5 7 ; 2 4 6 8];

figure 
subplot(2,1,1)
hold on
title('staccato')
plot(t,All_data(grp(1,d)).suj(suj).mu(ordre_mu(mu)).data_selec_10')

subplot(2,1,2)
hold on
title('tenuto')
plot(t,All_data(grp(2,d)).suj(suj).mu(ordre_mu(mu)).data_selec_10')

%% visualisation variabilité 
% body-implication

mu = 4;
ordre_mu = [9 8 2 1 3 4 7 5 6];
F=2100/5;
t = 1/F -630/(2*F) : 1/F : 630/(2*F);

m_b = mean( mean( MU(mu).boxplot_data(:,[1 2 3 4])));
m_ms = mean( mean( MU(mu).boxplot_data(:,[5 6 7 8])));

[a b] = min(abs(MU(mu).boxplot_data(:, [1 2 3 4]) - m_b) + abs(MU(mu).boxplot_data(:, [5 6 7 8]) - m_ms));
[c d] = min(a);
suj = b(d);
grp = [1 2 3 4 ; 5 6 7 8];

figure 
subplot(2,1,1)
hold on
title('all-body')
plot(t,All_data(grp(1,d)).suj(suj).mu(ordre_mu(mu)).data_selec_10')

subplot(2,1,2)
hold on
title('upper-limb')
plot(t,All_data(grp(2,d)).suj(suj).mu(ordre_mu(mu)).data_selec_10')

%% comparaison CV_inter

sujets = [];
Y = [];
for i_sub = 1:8
   Y = [Y ,[All_data(i_sub).suj_mean(ordre_mu).VR_inter]'];
end
%%
suj_set = [1:12]';
Surface = Surface'
iterations = 2000 ; %1000;
A = [zeros(12,1); zeros(12,1); ones(12,1); ones(12,1); ...
    zeros(12,1); zeros(12,1); ones(12,1); ones(12,1)]; % frappé (0) / pressé (1)

B = [zeros(12,1); zeros(12,1); zeros(12,1); zeros(12,1); ...
    ones(12,1); ones(12,1); ones(12,1); ones(12,1)]; % Bassin (0) / MS (1)

C = [zeros(12,1); ones(12,1); zeros(12,1); ones(12,1); ...
    zeros(12,1); ones(12,1); zeros(12,1); ones(12,1)]; % Sta (0) / Tenue (1)


Subj = [suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set ;suj_set];

F = spm1d.stats.nonparam.anova3rm(Surface(:), A, B, C,Subj);
Fi =F.inference(0.05,'iteration',iterations);
disp_summ(Fi)

%%
i_sub=0;
figure(1)
for i = 1 : 4
    i_sub = i_sub+1;
    subplot(2,2,i_sub)
    radarplot(Y(:,[i i+4])',ordre_mu_name, {'r' 'b'})
end

i_sub=0;
figure(2)
for i = [1 2 5 6]
    i_sub = i_sub+1;
    subplot(2,2,i_sub)
    radarplot(Y(:,[i i+2])',ordre_mu_name, {'r' 'b'})
end

i_sub=0;
figure(3)
for i = [1 3 5 7]
    i_sub = i_sub+1;
    subplot(2,2,i_sub)
    radarplot(Y(:,[i i+1])',ordre_mu_name, {'r' 'b'})
end


%% spider plot
close all
fig = figure;

grp1 = [1 2 5 6];
grp2 = [3 4 7 8];

R1 = mean(Y(:,grp1)') + std(Y(:,grp1)');
n1 = size(R1,2);
R1=[R1 R1(:,1)];

[Theta,M]=meshgrid(2*pi/n1*[0:n1]+pi/n1,ones(1,size(R1,1)));
X1=R1.*sin(Theta);
Y1=R1.*cos(Theta);

R2 = mean(Y(:,grp1)') - std(Y(:,grp1)');
R2=[R2 R2(:,1)];
[Theta,M]=meshgrid(2*pi/n1*[0:n1]+pi/n1,ones(1,size(R2,1)));
X2=R2.*sin(Theta);
Y2=R2.*cos(Theta);

hold on;
%     F1=fill(X1,Y1,cell2mat(grp_color(grp)),'LineStyle','none');
%     F2=fill(X2,Y2,cell2mat({'w'}),'LineStyle','none');
%     set(F1,'FaceAlpha',0.1)
%     set(F2,'FaceAlpha',1)
pa = patch([X1 X2], [Y1 Y2], 'k');
set(pa,'FaceAlpha',0.15, 'EdgeColor','None' )
set(get(get(pa,'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');   

%   %   %   %   %

R1 = mean(Y(:,grp2)') + std(Y(:,grp2)');

R1=[R1 R1(:,1)];
[Theta,M]=meshgrid(2*pi/n1*[0:n1]+pi/n1,ones(1,size(R1,1)));
X1=R1.*sin(Theta);
Y1=R1.*cos(Theta);

R2 = mean(Y(:,grp2)') - std(Y(:,grp2)');
R2=[R2 R2(:,1)];
[Theta,M]=meshgrid(2*pi/n1*[0:n1]+pi/n1,ones(1,size(R2,1)));
X2=R2.*sin(Theta);
Y2=R2.*cos(Theta);

hold on;
%     F1=fill(X1,Y1,cell2mat(grp_color(grp)),'LineStyle','none');
%     F2=fill(X2,Y2,cell2mat({'w'}),'LineStyle','none');
%     set(F1,'FaceAlpha',0.1)
%     set(F2,'FaceAlpha',1)
pa = patch([X1 X2], [Y1 Y2], 'k');
set(pa,'FaceAlpha',0.15, 'EdgeColor','None' )
set(get(get(pa,'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off'); 

radarplot([mean(Y(:,grp1)') ;mean(Y(:,grp2)')],ordre_mu_name)
legend({'Y' 'Y'})

set(gcf, 'Position', get(0, 'Screensize'))
%saveas(figure(1), '\\10.89.24.15\j\Valentin\Article EMG\figures\spider_exemple', 'jpeg')
%%

grp1 = [1 3 5 7];
grp2 = [2 4 6 8];

[h,p] = ttest2(Y(:, grp1)', Y(:, grp2)');
%% Cinematique de la touche


for grp = 1 : 8
   [S, I] = sort(mean(GRP(grp).boxplot_data));
%    [S, I] = sort(GRP(grp).boxplot_data_VR_inter);
    for mu = 1 : 9
        ind_mu_inter(grp,mu) = find(I == mu);
    end
%     hold on
%     plot(GRP(grp).boxplot_data_VR_inter)
    subplot(2,4,grp)
    boxplot(GRP(grp).boxplot_data(:,I), ordre_mu_name(I))
    ylim([0 1])
end

%%
for grp = 1 : 8
    b(grp,:) = GRP(grp).boxplot_data_VR_inter;
end



body = [1 1 1 1 0 0 0 0];
touch = [1 1 0 0 1 1 0 0];
arti = [1 0 1 0 1 0 1 0];

mean(b(arti==1,:))
std(b(arti==1,:))

%% spider plot
close all
fig = figure;

grp1 = [1 3 5 7];
grp2 = [2 4 6 8];
Y1_data = [];
Y2_data = [];
for grp = grp1
   Y1_data = [Y1_data ;GRP(grp).boxplot_data];
end

for grp = grp2
    Y2_data = [Y2_data ;GRP(grp).boxplot_data]; 
end
R1 = mean(Y1_data) + std(Y1_data);
n1 = size(R1,2);
R1=[R1 R1(:,1)];
[Theta,M]=meshgrid(2*pi/n1*[0:n1]+pi/n1,ones(1,size(R1,1)));
X1=R1.*sin(Theta);
Y1=R1.*cos(Theta);

R2 = mean(Y1_data) - std(Y1_data);
R2=[R2 R2(:,1)];
[Theta,M]=meshgrid(2*pi/n1*[0:n1]+pi/n1,ones(1,size(R2,1)));
X2=R2.*sin(Theta);
Y2=R2.*cos(Theta);

hold on;
%     F1=fill(X1,Y1,cell2mat(grp_color(grp)),'LineStyle','none');
%     F2=fill(X2,Y2,cell2mat({'w'}),'LineStyle','none');
%     set(F1,'FaceAlpha',0.1)
%     set(F2,'FaceAlpha',1)
pa = patch([X1 X2], [Y1 Y2], 'k');
set(pa,'FaceAlpha',0.15, 'EdgeColor','None' )
set(get(get(pa,'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');   

%   %   %   %   %

R1 = mean(Y2_data) + std(Y2_data);

R1=[R1 R1(:,1)];
[Theta,M]=meshgrid(2*pi/n1*[0:n1]+pi/n1,ones(1,size(R1,1)));
X1=R1.*sin(Theta);
Y1=R1.*cos(Theta);

R2 = mean(Y2_data) - std(Y2_data);
R2=[R2 R2(:,1)];
[Theta,M]=meshgrid(2*pi/n1*[0:n1]+pi/n1,ones(1,size(R2,1)));
X2=R2.*sin(Theta);
Y2=R2.*cos(Theta);

hold on;
%     F1=fill(X1,Y1,cell2mat(grp_color(grp)),'LineStyle','none');
%     F2=fill(X2,Y2,cell2mat({'w'}),'LineStyle','none');
%     set(F1,'FaceAlpha',0.1)
%     set(F2,'FaceAlpha',1)
pa = patch([X1 X2], [Y1 Y2], 'k');
set(pa,'FaceAlpha',0.15, 'EdgeColor','None' )
set(get(get(pa,'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off'); 

radarplot([mean(Y1_data) ;mean(Y2_data)],ordre_mu_name)
legend({'Y' 'Y'})


set(gcf, 'Position', get(0, 'Screensize'))
saveas(figure(1), '\\10.89.24.15\j\Valentin\Article EMG\figures\spider_surface_ssregression2', 'jpeg')
%%

F=150;
t = 1/F -227/(2*F) : 1/F : 227/(2*F);
Col = [0.87 0.49 0 ;1 0.84 0; 0 0 0; 0.49 0.49 0.49; 0.87 0.49 0 ;1 0.84 0; 0 0 0; 0.49 0.49 0.49];
% for grp = 1:8
%     hold on
%     err_c1 = spm1d.plot.plot_errorcloud(t, mean(Cine_touche(grp).grp), std(Cine_touche(grp).grp), 'facecolor', Col(grp,:), 'facealpha',0.20) ;
%     err_c1.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
% end
sty = {'-' '-' '-' '-' '--' '--' '--' '--'};
for grp = 1:8
    hold on
    plot(t, mean(Cine_touche(grp).grp),'Color', Col(grp,:),'Linestyle', sty{grp}, 'linewidth',2) ;
end

hold on
line([0 0], [-9 1])

%%
grp = 1
for mu = 1 : 9
    figure(mu)
    plot(All_data(grp).suj_mean(mu).mu') 
end

%% spider plot
close all

for grp = 1:8
    Y1_data = [];
    Y1_data = [Y1_data ;GRP(grp).boxplot_data];
    for suj = 1 : 12
        R1 = Y1_data(suj,:);
        n1 = size(R1,2);
        R1=[R1 R1(:,1)];
        [Theta,M]=meshgrid(2*pi/n1*[0:n1]+pi/n1,ones(1,size(R1,1)));
        X1_m=R1.*sin(Theta);
        Y1_m=R1.*cos(Theta);
    
    Surface(grp,suj) = polyarea(X1_m, Y1_m);
    end
end

%%
for grp = 1:8
    Y1_data = [];
    Y1_data = GRP(grp).boxplot_data_VR_inter;

    R1 = Y1_data;
    n1 = size(R1,2);
    R1=[R1 R1(:,1)];
    [Theta,M]=meshgrid(2*pi/n1*[0:n1]+pi/n1,ones(1,size(R1,1)));
    X1=R1.*sin(Theta);
    Y1=R1.*cos(Theta);

    Surface_inter(grp) = polyarea(X1, Y1);

end

[q,I]=sort(mean(Surface'));
figure(1)
hold on 
for grp = 1:8
    plot(grp, Surface_inter(I(grp)), 'o')
end
boxplot(Surface(I,:)',condition_name(I))
ylim([0 3])
% 
% hold on
% 
% x = 1:8;
% plot(x,sort(mean(Surface')), 'o')
% b1 = x'\sort(mean(Surface')-min(mean(Surface')))';
% plot(x, x*b1+min(mean(Surface')));
% 
% 
% plot(Surface_inter(I), 'p')
% b2 = x'\(Surface_inter(I)-min(Surface_inter(I)))';
% plot(x, x*b2+min(Surface_inter(I)));
% 
% corrcoef(sort(mean(Surface')), Surface_inter(I))
