clear all;
close all;
clc;


%addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))

pathname='C:\Users\p1218107\Documents\Data_Piano';

ind_suj = [2:7 9:12]; 
ind_mu = [1 4 5 6 9];

%% Ouverture des datas 

cd(pathname)

load([pathname,'\data_1.5s\Post_ttt_semiauto.mat'])

t0= round(length(Post_ttt(1).suj_mean(1).data)/2);

F=2100/4;

t = round((1-t0)/F: 1/F : (t0-1)/F,4);
%%

parametre = [];


for grp = 1 : 8 
    i=0;
    for mu = ind_mu
        
        if mu == 1 
            i = i+1;
            Data = (Post_ttt(grp).suj_mean(1).data + Post_ttt(grp).suj_mean(2).data)/2;
            parametre(grp).EMG_t0(:,i) = Data(ind_suj,t0);
            [parametre(grp).EMG_max(:,i), parametre(grp).tpic(:,i)] = max(Data(ind_suj,:), [],2);
        else
            i = i+1;
            parametre(grp).EMG_t0(:,i) = Post_ttt(grp).suj_mean(mu).data(ind_suj,t0);
            [parametre(grp).EMG_max(:,i), parametre(grp).tpic(:,i)] = max(Post_ttt(grp).suj_mean(mu).data(ind_suj,:), [],2);
        end
    end
    parametre(grp).tpic(parametre(grp).tpic>600)=NaN;
    parametre(grp).tpic(parametre(grp).tpic<200)=NaN;
end

%save([pathname, '\data_1.5s\', 'parametre'], 'parametre')

%% 
 Grp_names = {'Membre supérieur pressé staccato',...
  'Membre supérieur frappé staccato',...
  'Bassin pressé staccato',...
  'Bassin frappé staccato',...
  'Membre supérieur pressé tenue',...
  'Membre supérieur frappé tenue',...
  'Bassin pressé tenue'...
  'Bassin frappé tenue'};

for grp = 1 : 8
    subplot(4,2,grp)
    [B,I] = sort( nanmedian( parametre(grp).tpic));
    boxplot(parametre(grp).tpic(:, I) )
    hold on 
    line(xlim, [t0 t0], 'linestyle', '--', 'Color', [1 0 0] )
    title(Grp_names{grp})
    Mu_Names = {'Fléchisseurs';'Biceps';'Triceps';'DeltAnt';'TrapSup'};
    set(gca, 'xticklabel', Mu_Names(I))
end


% Data = [];
% for grp = 1 : 8
%     Data = [Data, parametre(grp).EMG_t0, parametre(grp).EMG_max, parametre(grp).tpic];
% end
% 
% mapcaplot(Data)
% 
% %%
% [coeff,score,latent,tsquared,explained] = pca(Data);
% 
% scatter3(score(:,1),score(:,2),score(:,3))
% axis equal
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')
% zlabel('3rd Principal Component')

