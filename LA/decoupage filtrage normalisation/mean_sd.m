%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        Traitement des EMGs                              %
%                          Par statistique                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation

clear all;
close all;
clc;

%addpath(genpath('C:\Users\Deg\Documents\MATLAB\ezc3d'))
addpath(genpath('C:\Users\p1218107\Documents\MATLAB\spm1dmatlab-master'))
sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};

pathname='C:\Users\p1218107\Documents\Data_Piano\data_1.5s\'; 





%% Sélection des Datas à traiter

for n = 1

    keepvars = {'n','nbframes' , 'pathname' , 'sujet'};
    
    close all

    
    clearvars('-except', keepvars{:});
    
    PathName_EMG=[pathname,sujet{n},'_EMG_GRP.mat'];

    load(PathName_EMG)

    Names={'Groupe 1','Groupe 2','Groupe 3','Groupe 4','Groupe 5','Groupe 6','Groupe 7','Groupe 8'};
    %% TTT GRPS

    tic
    for grp = 4 %1:4 

        %%
            for mu = 10 % 1 : 10
                

                h = figure(grp);
                %set(gcf, 'Position', get(0, 'Screensize'));
                name=[sujet{n},' ',GRP(grp).Name];
                suptitle(name)
                %subplot(5,2,mu)
                for cycl = 1 : size(GRP(grp).EMGred(mu).data_norm,1)
                    hold on
                    plot(GRP(grp).EMGred(mu).data_norm(cycl,:),'k')
                end
                meanA = [];
                meanA = mean(GRP(grp).EMGred(mu).data_norm,1);
                
                sd = std(GRP(grp).EMGred(mu).data_norm);

                plot(meanA+3*sd,'r')
                %plot(meanA-3*sd,'b')
                
                
                B=(meanA+3*sd);
                
                GRP(grp).EMGred(mu).data_ttt2 = GRP(grp).EMGred(mu).data_norm ; 
                %%
                
                for realign = 1:4
                    for ttt = 1:10
                    A=[];
                    o=0;
                        for cycl = 1 : size(GRP(grp).EMGred(mu).data_ttt2,1)

                            if sum(GRP(grp).EMGred(mu).data_ttt2(cycl,[365 425])>B([365 425]))/length(B)==0%<0.1
                                o=o+1;
                                A(o,:)=GRP(grp).EMGred(mu).data_norm(cycl,:);
%                                 hold on
%                                 color={'g','m','c','r'};
%                                 plot(GRP(grp).EMGred(mu).data(cycl,:),color{ttt})
                            end
                        end

                        meanA= mean(A,1);
                        sd = std(A);
                        B=(meanA+3*sd);

                    end

                    for cycl = 1 : size(A,1)

                        D = finddelay(A(cycl,100:end), meanA);
                            
                        if D < 10
                            A(cycl,:) =  A(cycl,:); 
                        end
                        if D<99 && D > 10 
                            A(cycl,:) = [A(cycl,100-D:end), meanA(1,1:99-D)];
                        end
                        if D>=100
                            A(cycl,:) = [meanA(1,1:D-99), A(cycl,1:end+99-D)];
                        end

                    end

                    meanA= mean(A,1);
                    sd = std(A);
                    B=(meanA+3*sd);

                    GRP(grp).EMGred(mu).data_ttt2=A;
                end
                
                GRP(grp).EMGred(mu).data_norm = GRP(grp).EMGred(mu).data_ttt2;

                
                subplot(5,2,mu)
                for cycl = 1 : size(GRP(grp).EMGred(mu).data_ttt2,1)
                    hold on
                    plot(GRP(grp).EMGred(mu).data_ttt2(cycl,:),'y')
                end
                
                plot(meanA,'r','linewidth',1.5)
                title(GRP(1).EMGred(mu).labels)

            end
            
            %saveas(h,['C:\Users\p1218107\Documents\Data_Piano\graphiques\ttt_stat\',sujet{n},'_grp',num2str(grp)],'jpg')
            
    end
    %save(PathName_EMG,'GRP')
    toc
end


