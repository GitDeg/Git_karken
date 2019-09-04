clear all;
close all;
clc;

addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};

%%
direction='C:\Users\Deg\Desktop\Thèse\DATA_Piano';

REPS={};

for n=1:length(sujet)
    spe_dir=[direction,'\',sujet{n},'\LA','\all_EMG_',sujet{n},'_cycle.mat'];
    
    load(spe_dir);
    load([direction,'\',sujet{n},'\LA','\all_Fo_',sujet{n},'.mat'])
    
    for i = 1:length(EMG_LA)
        REPS{i,n*4-3}=EMG_LA(i).reps;
        REPS{i,n*4-2}=EMG_LA(i).Filename;
        REPS{i,n*4-1}=EMG_LA(i).Rate;
        REPS{i,n*4}=length(Fo_LA(i).RawData);
    end
    
    clear('EMG_LA','Fo_LA')
end

%%

save('all_frames.mat','REPS');