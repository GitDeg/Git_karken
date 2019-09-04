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

for n = 6
   
    

    keepvars = {'n','nbframes' , 'pathname' , 'sujet'};
    
    close all

    
    clearvars('-except', keepvars{:});
    
    PathName_EMG=[pathname,sujet{n},'_EMG_GRP.mat'];

    load(PathName_EMG)

    Names={'Groupe 1','Groupe 2','Groupe 3','Groupe 4','Groupe 5','Groupe 6','Groupe 7','Groupe 8'};
    %% TTT GRPS

    tic
    for grp = 5:8

        %%
            for mu = 1:10
                
                close all
                
                GRP(grp).EMGred(mu).data_ttt3 = GRP(grp).EMGred(mu).data ;
                %%
                ttt = 0;
                while 1
                    
                    ttt = ttt + 1 ;
                    
                    h = figure(ttt);
                    set(gcf, 'Position', get(0, 'Screensize'))

                    %set(gcf, 'Position', get(0, 'Screensize'));
                    name=[sujet{n},' ',GRP(grp).Name];
                    suptitle(name)
                    %subplot(5,2,mu)
                    for cycl = 1 : size(GRP(grp).EMGred(mu).data_ttt3,1)
                        hold on
                        plot(GRP(grp).EMGred(mu).data_ttt3(cycl,:),'k')
                    end
                    meanA = [];
                    meanA = mean(GRP(grp).EMGred(mu).data_ttt3,1);

                    sd = std(GRP(grp).EMGred(mu).data_ttt3);

                    plot(meanA+3*sd,'r')
                    plot(meanA, 'g', 'linewidth', 2)
                    %plot(meanA-3*sd,'b')


                    B=(meanA+3*sd);
                    dcm_obj = datacursormode(h);
                    set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on')
                    prompt = ['TTT? Y/N: '];
                    str = input(prompt,'s');
                    
                    if str=='N'
                        break
                    else 
                        c_info = getCursorInfo(dcm_obj);
                        
                        bind = round(c_info.Position(1));
                    end
                    %%
                    %[bmax,bind] = max(B);
                    
                    var_evo = [];
                    for cycl = 1 : size(GRP(grp).EMGred(mu).data_ttt3,1)
                        A = 1 : size(GRP(grp).EMGred(mu).data_ttt3,1);
                        var_evo(cycl) = std(GRP(grp).EMGred(mu).data_ttt3(A(A~=cycl),bind));
                    end
                    
%                     J = figure(50+ttt);
%                     set(gcf, 'Position', get(0, 'Screensize'))
% 
%                     plot(var_evo)
                    
%                     dcm_obj = datacursormode(J);
%                     set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on')
%                     
%                     c_info = [] ; 
                    
%                     
%                     waitforbuttonpress
%                     waitforbuttonpress
%                     c_info = getCursorInfo(dcm_obj);
%                     
% 
%                     seuil = c_info.Position(2);

                    [a,b] = min(var_evo);
                    seuil = a;
                    x_var = [];
                    x_var = b;
%                     x_var = find(var_evo<seuil);

                    GRP(grp).EMGred(mu).data_ttt3(x_var,:)=[];
                    figure(ttt)
                    for cycl = 1 : size(GRP(grp).EMGred(mu).data_ttt3,1)
                        hold on
                        plot(GRP(grp).EMGred(mu).data_ttt3(cycl,:),'y')
                    end
                    
                end

            end
            
            %saveas(h,['C:\Users\p1218107\Documents\Data_Piano\graphiques\ttt_stat\',sujet{n},'_grp',num2str(grp)],'jpg')
            
    end
    save(PathName_EMG,'GRP')
    toc
end


