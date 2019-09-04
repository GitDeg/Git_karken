%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Chargement des datas EMG et Fo et trie en 4 groupes           %                                       %
%                              Protocole LA                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation

clear all;
close all;
clc;

addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
pathname='C:\Users\Deg\Desktop\Thèse\DATA_Piano';
sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};

%% Sélection des Datas à traiter


for n=1:12
%     cd(pathname)
%     [FileName_EMG,PathName] = uigetfile('*.mat',['Selection Data all_EMG_',num2str(n),'cycle','.mat']);
%     cd(PathName)
    load([pathname,'\',sujet{n},'\LA\all_EMG_',sujet{n},'_cycle.mat'])

    load([pathname,'\',sujet{n},'\LA\all_Fo_',sujet{n},'.mat'])
    %EMG_LAf=[];

%% Tri en 8 groupes:
%
% Bas = Tout le corps 
% Sup = Membre supérieur 
% Pre = Toucher pressé 
% Fra = Toucher frappé 
% Ten = Note tenue 
% Sta = Note stacato 
%
% Groupe 1: Membre supérieur pressé staccato -> Grp1
% Groupe 2: Membre supérieur frappé staccato -> Grp2
% Groupe 3: Bassin pressé staccato -> Grp3
% Groupe 4: Bassin frappé staccato -> Grp4
% Groupe 5: Membre supérieur pressé tenue ->Grp5
% Groupe 6: Membre supérieur frappé tenue ->Grp6
% Groupe 7: Bassin pressé tenue ->Grp7
% Groupe 8: Bassin frappé tenue ->Grp8

    Grps={'Grp1','Grp2','Grp3','Grp4','Grp5','Grp6','Grp7','Grp8'};

    Grp1.EMG= [];
    Grp1.EMG.RawData= [];
    Grp1.EMG.Labels= {};
    Grp1.EMG.Rate= [];
    Grp1.EMG.Filename= [];
    Grp1.EMG.reps= [];
    Grp1.Fo.RawData= [];
    Grp1.Fo.Labels= [];


    Grp2.EMG= [];
    Grp2.EMG.RawData= [];
    Grp2.EMG.Labels= {};
    Grp2.EMG.Rate= [];
    Grp2.EMG.Filename= [];
    Grp2.EMG.reps= [];
    Grp2.Fo.RawData= [];
    Grp2.Fo.Labels= [];

    Grp3.EMG= [];
    Grp3.EMG.RawData= [];
    Grp3.EMG.Labels= {};
    Grp3.EMG.Rate= [];
    Grp3.EMG.Filename= [];
    Grp3.EMG.reps= [];
    Grp3.Fo.RawData= [];   
    Grp3.Fo.Labels= [];

    Grp4.EMG= [];
    Grp4.EMG.RawData= [];
    Grp4.EMG.Labels= {};
    Grp4.EMG.Rate= [];
    Grp4.EMG.Filename= [];
    Grp4.EMG.reps= [];
    Grp4.Fo.RawData= [];
    Grp4.Fo.Labels= [];
    
    
    Grp5.EMG= [];
    Grp5.EMG.RawData= [];
    Grp5.EMG.Labels= {};
    Grp5.EMG.Rate= [];
    Grp5.EMG.Filename= [];
    Grp5.EMG.reps= [];
    Grp5.Fo.RawData= [];
    Grp5.Fo.Labels= [];


    Grp6.EMG= [];
    Grp6.EMG.RawData= [];
    Grp6.EMG.Labels= {};
    Grp6.EMG.Rate= [];
    Grp6.EMG.Filename= [];
    Grp6.EMG.reps= [];
    Grp6.Fo.RawData= [];
    Grp6.Fo.Labels= [];

    Grp7.EMG= [];
    Grp7.EMG.RawData= [];
    Grp7.EMG.Labels= {};
    Grp7.EMG.Rate= [];
    Grp7.EMG.Filename= [];
    Grp7.EMG.reps= [];
    Grp7.Fo.RawData= [];   
    Grp7.Fo.Labels= [];

    Grp8.EMG= [];
    Grp8.EMG.RawData= [];
    Grp8.EMG.Labels= {};
    Grp8.EMG.Rate= [];
    Grp8.EMG.Filename= [];
    Grp8.EMG.reps= [];
    Grp8.Fo.RawData= [];
    Grp8.Fo.Labels= [];


%% Boucle de tri + suppression des données inutiles

    SPSi= 0;
    SFSi= 0;
    BPSi= 0;
    BFSi= 0;
    
    SPTi= 0;
    SFTi= 0;
    BPTi= 0;
    BFTi= 0;

%%
    for i=1:length(EMG_LA)
   
        if all(ismember('Ten',EMG_LA(i).Filename))==1 && all(ismember(['Bas','Pre'],EMG_LA(i).Filename))==1
            if BPTi==0
                BPTi=BPTi+1;
                Grp7.EMG(BPTi)=EMG_LA(i);
                Grp7.Fo(BPTi).RawData= Fo_LA(i).RawData;
                Grp7.Fo(BPTi).Labels= Fo_LA(i).Labels;
            elseif ismember(EMG_LA(i).Filename(10),Grp7.EMG(BPTi).Filename(10))==1 %si le nom existe déjà, il remplace les valeurs enregistrés.
                Grp7.EMG(BPTi)=EMG_LA(i);
                Grp7.Fo(BPTi).RawData= Fo_LA(i).RawData;
                Grp7.Fo(BPTi).Labels= Fo_LA(i).Labels;
            else
                BPTi=BPTi+1;
                Grp7.EMG(BPTi)=EMG_LA(i); 
                Grp7.Fo(BPTi).RawData= Fo_LA(i).RawData;
                Grp7.Fo(BPTi).Labels= Fo_LA(i).Labels;
            end
        end
        
        if all(ismember('Ten',EMG_LA(i).Filename))==1 && all(ismember(['Bas','Fra'],EMG_LA(i).Filename))==1
            if BFTi==0
                BFTi=BFTi+1;
                Grp8.EMG(BFTi)=EMG_LA(i);
                Grp8.Fo(BFTi).RawData= Fo_LA(i).RawData;
                Grp8.Fo(BFTi).Labels= Fo_LA(i).Labels;
            elseif ismember(EMG_LA(i).Filename(10),Grp8.EMG(BFTi).Filename(10))==1
                Grp8.EMG(BFTi)=EMG_LA(i);
                Grp8.Fo(BFTi).RawData= Fo_LA(i).RawData;
                Grp8.Fo(BFTi).Labels= Fo_LA(i).Labels;
            else
                BFTi=BFTi+1;
                Grp8.EMG(BFTi)=EMG_LA(i); 
                Grp8.Fo(BFTi).RawData= Fo_LA(i).RawData;
                Grp8.Fo(BFTi).Labels= Fo_LA(i).Labels;
            end
        end
        
        if all(ismember('Ten',EMG_LA(i).Filename))==1 && all(ismember(['Sup','Fra'],EMG_LA(i).Filename))==1
            if SFTi==0
                SFTi=SFTi+1;
                Grp6.EMG(SFTi)=EMG_LA(i);
                Grp6.Fo(SFTi).RawData= Fo_LA(i).RawData;
                Grp6.Fo(SFTi).Labels= Fo_LA(i).Labels;
            elseif ismember(EMG_LA(i).Filename(10),Grp6.EMG(SFTi).Filename(10))==1
                Grp6.EMG(SFTi)=EMG_LA(i);
                Grp6.Fo(SFTi).RawData= Fo_LA(i).RawData;
                Grp6.Fo(SFTi).Labels= Fo_LA(i).Labels;
            else
                SFTi=SFTi+1;
                Grp6.EMG(SFTi)=EMG_LA(i);  
                Grp6.Fo(SFTi).RawData= Fo_LA(i).RawData;
                Grp6.Fo(SFTi).Labels= Fo_LA(i).Labels;
            end
        end
        
        if all(ismember('Ten',EMG_LA(i).Filename))==1 && all(ismember(['Sup','Pre'],EMG_LA(i).Filename))==1
            if SPTi==0
                SPTi=SPTi+1;
                Grp5.EMG(SPTi)=EMG_LA(i);
                Grp5.Fo(SPTi).RawData= Fo_LA(i).RawData;
                Grp5.Fo(SPTi).Labels= Fo_LA(i).Labels;
            elseif ismember(EMG_LA(i).Filename(10),Grp5.EMG(SPTi).Filename(10))==1
                Grp5.EMG(SPTi)=EMG_LA(i);
                Grp5.Fo(SPTi).RawData= Fo_LA(i).RawData;
                Grp5.Fo(SPTi).Labels= Fo_LA(i).Labels;
            else
                SPTi=SPTi+1;
                Grp5.EMG(SPTi)=EMG_LA(i); 
                Grp5.Fo(SPTi).RawData= Fo_LA(i).RawData;
                Grp5.Fo(SPTi).Labels= Fo_LA(i).Labels;
            end
        end
                
        if all(ismember('Ten',EMG_LA(i).Filename))==0 && all(ismember(['Bas','Pre'],EMG_LA(i).Filename))==1
            if BPSi==0
                BPSi=BPSi+1;
                Grp3.EMG(BPSi)=EMG_LA(i);
                Grp3.Fo(BPSi).RawData= Fo_LA(i).RawData;
                Grp3.Fo(BPSi).Labels= Fo_LA(i).Labels;
            elseif ismember(EMG_LA(i).Filename(10),Grp3.EMG(BPSi).Filename(10))==1 %si le nom existe déjà, il remplace les valeurs enregistrés.
                Grp3.EMG(BPSi)=EMG_LA(i);
                Grp3.Fo(BPSi).RawData= Fo_LA(i).RawData;
                Grp3.Fo(BPSi).Labels= Fo_LA(i).Labels;
            else
                BPSi=BPSi+1;
                Grp3.EMG(BPSi)=EMG_LA(i); 
                Grp3.Fo(BPSi).RawData= Fo_LA(i).RawData;
                Grp3.Fo(BPSi).Labels= Fo_LA(i).Labels;
            end
        end
    
        if all(ismember('Ten',EMG_LA(i).Filename))==0 && all(ismember(['Bas','Fra'],EMG_LA(i).Filename))==1
            if BFSi==0
                BFSi=BFSi+1;
                Grp4.EMG(BFSi)=EMG_LA(i);
                Grp4.Fo(BFSi).RawData= Fo_LA(i).RawData;
                Grp4.Fo(BFSi).Labels= Fo_LA(i).Labels;
            elseif ismember(EMG_LA(i).Filename(10),Grp4.EMG(BFSi).Filename(10))==1
                Grp4.EMG(BFSi)=EMG_LA(i);
                Grp4.Fo(BFSi).RawData= Fo_LA(i).RawData;
                Grp4.Fo(BFSi).Labels= Fo_LA(i).Labels;
            else
                BFSi=BFSi+1;
                Grp4.EMG(BFSi)=EMG_LA(i); 
                Grp4.Fo(BFSi).RawData= Fo_LA(i).RawData;
                Grp4.Fo(BFSi).Labels= Fo_LA(i).Labels;
            end
        end
    
        if all(ismember('Ten',EMG_LA(i).Filename))==0 && all(ismember(['Sup','Fra'],EMG_LA(i).Filename))==1
            if SFSi==0
                SFSi=SFSi+1;
                Grp2.EMG(SFSi)=EMG_LA(i);
                Grp2.Fo(SFSi).RawData= Fo_LA(i).RawData;
                Grp2.Fo(SFSi).Labels= Fo_LA(i).Labels;
            elseif ismember(EMG_LA(i).Filename(10),Grp2.EMG(SFSi).Filename(10))==1
                Grp2.EMG(SFSi)=EMG_LA(i);
                Grp2.Fo(SFSi).RawData= Fo_LA(i).RawData;
                Grp2.Fo(SFSi).Labels= Fo_LA(i).Labels;
            else
                SFSi=SFSi+1;
                Grp2.EMG(SFSi)=EMG_LA(i);  
                Grp2.Fo(SFSi).RawData= Fo_LA(i).RawData;
                Grp2.Fo(SFSi).Labels= Fo_LA(i).Labels;
            end
        end
    
        if all(ismember('Ten',EMG_LA(i).Filename))==0 && all(ismember(['Sup','Pre'],EMG_LA(i).Filename))==1
            if SPSi==0
                SPSi=SPSi+1;
                Grp1.EMG(SPSi)=EMG_LA(i);
                Grp1.Fo(SPSi).RawData= Fo_LA(i).RawData;
                Grp1.Fo(SPSi).Labels= Fo_LA(i).Labels;
            elseif ismember(EMG_LA(i).Filename(10),Grp1.EMG(SPSi).Filename(10))==1
                Grp1.EMG(SPSi)=EMG_LA(i);
                Grp1.Fo(SPSi).RawData= Fo_LA(i).RawData;
                Grp1.Fo(SPSi).Labels= Fo_LA(i).Labels;
            else
                SPSi=SPSi+1;
                Grp1.EMG(SPSi)=EMG_LA(i); 
                Grp1.Fo(SPSi).RawData= Fo_LA(i).RawData;
                Grp1.Fo(SPSi).Labels= Fo_LA(i).Labels;
            end
        end
       
    end

    clear('SPSi', 'SFSi', 'BFSi');

%     save([pathname,'All_Grp1'],'Grp1')
%     save([pathname,'All_Grp2'],'Grp2')
%     save([pathname,'All_Grp3'],'Grp3')
%     save([pathname,'All_Grp4'],'Grp4')
    save([pathname,'\',sujet{n},'\LA\','All_Grp5'],'Grp5')
%    save([pathname,'\',sujet{n},'\LA\','All_Grp6'],'Grp6')
%    save([pathname,'\',sujet{n},'\LA\','All_Grp7'],'Grp7')
%    save([pathname,'\',sujet{n},'\LA\','All_Grp8'],'Grp8')
    
    clear('Grp1','Grp2','Grp3','Grp4','Grps')
    clear('Grp5','Grp6','Grp7','Grp8','Grps')
end
    %%
