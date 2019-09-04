clear all;
close all;
clc;


%addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB\spm1dmatlab-master'))

pathname='C:\Users\p1218107\Documents\Data_Piano';

%%
mu = 1;

% Y = [EMG_max(mu+1).grp1, EMG_contact(mu+1).grp1, tpic(mu+1).grp1; ...
%     EMG_max(mu+1).grp2, EMG_contact(mu+1).grp2, tpic(mu+1).grp2; ...
%     EMG_max(mu+1).grp3, EMG_contact(mu+1).grp3, tpic(mu+1).grp3; ...
%     EMG_max(mu+1).grp4, EMG_contact(mu+1).grp4, tpic(mu+1).grp4] ;

Y = [EMG_max(mu+1).grp1, EMG_max(mu+1).grp2, EMG_max(mu+1).grp3, EMG_max(mu+1).grp4] ;

clusterdata(Y,12);

idx = kmeans(Y , 12);

figure 
hold on

plot3(Y(idx ==1,1) , Y(idx == 1,2), Y(idx == 1,3), 'k*', 'MarkerSize', 5 );
plot3(Y(idx ==2,1) , Y(idx == 2,2), Y(idx == 2,3), 'r*', 'MarkerSize', 5 );
plot3(Y(idx ==3,1) , Y(idx == 3,2), Y(idx == 3,3), 'b*', 'MarkerSize', 5 );
plot3(Y(idx ==4,1) , Y(idx == 4,2), Y(idx == 4,3), 'g*', 'MarkerSize', 5 );

xlabel 'EMG max';
ylabel 'EMG contact';
zlabel 't pic';
grid on

%% ACP

suj = 1;
mu = 1;
l = length(para(mu).EMGmean(para(mu).suj==suj,:));
grp1= 1:l;
grp2= l + 1: 2*l;
grp3= 2*l + 1: 3*l;
grp4= 3*l + 1: 4*l;

Y = [reshape(para(mu).EMGmean(para(mu).suj==suj,:),[],1), reshape(para(mu).tpic(para(mu).suj==suj,:),[],1), reshape(para(mu).EMGp200(para(mu).suj==suj,:),[],1), ...
     reshape(para(mu).EMGmax(para(mu).suj==suj,:),[],1), reshape(para(mu).EMGmin(para(mu).suj==suj,:),[],1), reshape(para(mu).contact(para(mu).suj==suj,:),[],1),...
     reshape(para(mu).EMGm200(para(mu).suj==suj,:),[],1),  reshape(para(mu).EMGp50(para(mu).suj==suj,:),[],1),  reshape(para(mu).EMGm50(para(mu).suj==suj,:),[],1), ...
     reshape(para(mu).nbpic(para(mu).suj==suj,:),[],1)]; 

Z = zscore (Y);

[coeff,score,latent,tsquared,explained] = pca(Z);

scatter3(score(:,1),score(:,2),score(:,3))
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')



biplot(coeff(:,1:3),'varlabels',{'v_1','v_2','v_3','v_4','v_5','v_6','v_7','v_8','v_9','v_(10)'});
figure %% ,'scores',score(:,1:3)
hold on
grid on
plot3(score(grp1,1),score(grp1,2),score(grp1,3),'o','color','r')
plot3(score(grp2,1),score(grp2,2),score(grp2,3),'o','color','b')
plot3(score(grp3,1),score(grp3,2),score(grp3,3),'x','color','r')
plot3(score(grp4,1),score(grp4,2),score(grp4,3),'x','color','b')
mapcaplot(Z)

%%


idx = kmeans (score,4);

err(1,1) = sum(idx(grp1) == 1);
err(1,2) = sum(idx(grp1) == 2);
err(1,3) = sum(idx(grp1) == 3);
err(1,4) = sum(idx(grp1) == 4);

err(2,1) = sum(idx(grp2) == 1);
err(2,2) = sum(idx(grp2) == 2);
err(2,3) = sum(idx(grp2) == 3);
err(2,4) = sum(idx(grp2) == 4);

err(3,1) = sum(idx(grp3) == 1);
err(3,2) = sum(idx(grp3) == 2);
err(3,3) = sum(idx(grp3) == 3);
err(3,4) = sum(idx(grp3) == 4);

err(4,1) = sum(idx(grp4) == 1);
err(4,2) = sum(idx(grp4) == 2);
err(4,3) = sum(idx(grp4) == 3);
err(4,4) = sum(idx(grp4) == 4);




