clear all
close all
clc 

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
addpath(genpath('\\10.89.24.15\e\Librairies\S2M_Lib\Benjamin_dev\Plots'))
%%
% pathname = '\\10.89.24.15\j\Valentin\SAUT';
pathname = 'C:\Users\p1218107\Documents\Data_Piano\SAUT';
sujet = {'\001','\002','\003','\004','\005','\006','\007','\008','\009','\010','\011','\012'};
%%
suj = 6;

load([pathname sujet{suj} '\CINE_SAUT.mat'])
load([pathname sujet{suj} '\CINE_LA.mat'])

%% Repere piano
load([pathname sujet{suj} '\PIANO.mat'])
Gra_down = [];
Aig_down = [];
Aig_ant  = [];

for grp = [1 2 5 6]
    Gra_down = [Gra_down , Piano(grp).CINE(1).data];
    Aig_down = [Aig_down , Piano(grp).CINE(2).data];
    Aig_ant  = [Aig_ant  , Piano(grp).CINE(3).data];
end

Gra_down = nanmean(Gra_down')';
Aig_down = nanmean(Aig_down')';
Aig_ant  = nanmean(Aig_ant')';

%%
E4 = CINE_SAUT(1).RawData(2).Data(:,50:400);
E4(E4 == 0) = NaN;
E4 = nanmean(E4')';
A4 = CINE_SAUT(1).RawData(3).Data(:,50:400);
A4(A4 == 0) = NaN;
A4 = nanmean(A4')';
C5 = CINE_SAUT(1).RawData(4).Data(:,50:400);
C5(C5 == 0) = NaN;
C5 = nanmean(C5')';
C5_up = C5 + [0; 0; 100];
E5 = CINE_SAUT(1).RawData(5).Data(:,50:400);
E5(E5 == 0) = NaN;
E5 = nanmean(E5')';
C6 = CINE_SAUT(1).RawData(6).Data(:,50:400);
C6(C6 == 0) = NaN;
C6 = nanmean(C6')';

%%
X_axis = (Aig_down - Gra_down)/ norm(Aig_down - Gra_down) ;
v = (Aig_ant - Aig_down)/ norm(Aig_ant - Aig_down);
Z_axis = cross(X_axis,v)/ norm(cross(X_axis,v));
Y_axis = cross(Z_axis, X_axis);
P_G_P = [X_axis Y_axis Z_axis Gra_down; 0 0 0 1];

%% 
X_t_G = (C6(1:3) - E4(1:3));
v_z_G = (E4(1:3) - C5_up(1:3));
Y_t_G = cross(X_t_G, v_z_G);
Z_t_G = cross(X_t_G, Y_t_G);
X_t_G = X_t_G / norm(X_t_G);
Y_t_G = Y_t_G / norm(Y_t_G);
Z_t_G = Z_t_G / norm(Z_t_G);
P_G_T = [X_t_G Y_t_G Z_t_G C5; 0 0 0 1];

P_P_T = invR(P_G_P) * P_G_T;
P_G_T = P_G_P * P_P_T;

%%
figure 
hold on
axis equal
grid on

plotAxes(eye(4), 'length', 100); 
plotAxes(P_G_P, 'length', 100);
plotAxes(P_G_T, 'length', 100);
plot3d([Gra_down Aig_down Aig_ant], 'ko')
plot3d([E4 A4 C5 C5_up E5 C6], 'ro');
plot3d(E4, 'ro')
