%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Filtrage / Decoupage                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc 

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))

%%
pathname = '\\10.89.24.15\j\Valentin\SAUT';
sujet = {'\001','\002','\003','\004','\005','\006','\007','\008','\009','\010','\011','\012'};
%%
suj = 1;
load([pathname sujet{suj} '\CYCLES_SAUT.mat'])
%%
    mu = 10;
    Y = [];
    i = 0;
    figure
    for grp = [2 10]% 18]
        i=i+1;
        subplot(2,1,i)
        hold on
        plot(Data(grp).mu(mu).cycle_norm')
        plot(mean(Data(grp).mu(mu).cycle_norm), 'Linewidth', 2, 'Color', 'r')
        Y = [Y; Data(grp).mu(mu).cycle_norm]; 
    end
    
    A = [zeros(15,1); ones(15,1)];%; 2*ones(15,1)];
    
    figure
    F = spm1d.stats.nonparam.anova1(Y, A);
    Fi =F.inference(0.05,'iteration',1000);
    disp(Fi);
    plot(Fi);
    