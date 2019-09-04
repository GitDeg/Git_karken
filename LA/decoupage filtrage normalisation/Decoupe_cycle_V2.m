%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% redécoupe cycle data force
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation

clear all;
close all;
clc;

addpath(genpath('C:\Users\p1218107\Documents\MATLAB'))
sujet={'001','002','003','004','005','006','007','008','009','010','011','012'};
pathname = 'C:\Users\p1218107\Documents\Data_Piano\';
pathname_Cine = '\\10.89.24.15\e\Projet_Reconstructions\Piano_Reconstructions\script_matlab\';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                           Données Forces                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for suj = 1
    for grp = 1 : 8
        suff = num2str(grp);
        load([pathname sujet{suj} '\' sprintf('All_Grp%s',suff)]);
    end
%%
    GRPs(1).EMG = Grp1.EMG;
    GRPs(1).Fo = Grp1.Fo;
    
    GRPs(2).EMG = Grp2.EMG;
    GRPs(2).Fo = Grp2.Fo;
    
    GRPs(3).EMG = Grp3.EMG;
    GRPs(3).Fo = Grp3.Fo;
    
    GRPs(4).EMG = Grp4.EMG;
    GRPs(4).Fo = Grp4.Fo;
    
    GRPs(5).EMG = Grp5.EMG;
    GRPs(5).Fo = Grp5.Fo;
    
    GRPs(6).EMG = Grp6.EMG;
    GRPs(6).Fo = Grp6.Fo;
    
    GRPs(7).EMG = Grp7.EMG;
    GRPs(7).Fo = Grp7.Fo;
    
    GRPs(8).EMG = Grp8.EMG;
    GRPs(8).Fo = Grp8.Fo;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
    close all 
    
    for grp = 1 : 8
        % Force 1, premier serie d'essai 
        data = sum(GRPs(grp).Fo(1).RawData); 
        sdata = sort(data); 
        
        Force1 = 50 * (data - mean(sdata(1 : length(sdata)/4)));
        Force1p = sgolayfilt(sgolayfilt(Force1,1,49),1,49);



        %plot(Force1)

        Save = 'N';
        while Save == 'N'
            close all
            fig = figure('units','normalized','outerposition',[0 0 1 1]);
            plot(Force1p)
            hold on
            clear('locs', 'pks')

            prompt = 'Manuel(1) : choix du seuil et des pics / Automatique(2) : Choix du seuil large / Semi-Automatique(3) : Choix du seuil restreint ? [ 1 / 2 / 3 ]';
            M_A = input(prompt,'s');

            if M_A == '1'
                dcm_obj = datacursormode(fig);
                set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on')

                str='N';
                while str == 'N'
                    prompt ='Plot des points? [Y/N]';
                    str = input(prompt,'s');

                    if str ~= 'N'
                        c_info = getCursorInfo(dcm_obj);
                        hold on
                        pks = Force1p([c_info.DataIndex]);
                        locs = [c_info.DataIndex];
                        plot(locs, pks, 'o')
                        seuil = input('seuil ? ');
                        Force1_pics = Force1p;
                        Force1_pics(Force1_pics<seuil) = 0;
                        
                        delta = 0.1;
                        
                            for i = 1 : length(pks)
                                dif = Force1_pics(locs(i));
                                while dif > delta
                                    locs(i) = locs(i) - 2 ; 
                                    dif = Force1_pics(locs(i));
                                end
                            end
                            plot(locs , Force1(locs), 'ks')
                    end
                end
            end

            if M_A == '3'
                seuil = input('Valeur du Seuil ? ');
                line(xlim, [seuil seuil], 'Linestyle', '--', 'Color', [1 0 0])

                Force1_pics = Force1p;
                Force1_pics(Force1_pics<seuil) = 0; 

                [pks, locs] = findpeaks(Force1_pics);
                plot(locs, pks, 'o')

                delta = 0.1;

                for i = 1 : length(pks)

                    dif = Force1_pics(locs(i));
                    while dif > delta
                        locs(i) = locs(i) - 2 ; 
                        dif = Force1_pics(locs(i));
                    end
                end
                plot(locs , Force1(locs), 'ks')
            end
            
            if M_A == '2'
                seuil = input('Valeur du Seuil ? ');
                line(xlim, [seuil seuil], 'Linestyle', '--', 'Color', [1 0 0])

                Force1_pics = Force1p;
                Force1_pics(Force1_pics<seuil) = 0; 

                [pks, locs] = findpeaks(Force1_pics);
                plot(locs, pks, 'o')

                delta = 0.1;

                for i = 1 : length(pks)

                    dif = Force1p(locs(i));
                    while dif > delta
                        locs(i) = locs(i) - 2 ; 
                        dif = Force1p(locs(i));
                    end
                end
                plot(locs , Force1(locs), 'ks')
            end

            Save = input('Save ? [ Y / N ] ', 's');

            if Save ~= 'N'
                GRPs(grp).EMG(1).reps2 = locs(find([100 diff(locs)]>99));
            end
        end

    %%   

        % Force 2, deuxieme serie d'essai 
        
        data = sum(GRPs(grp).Fo(2).RawData); 
        sdata = sort(data); 
        
        Force2 = 50 * (data - mean(sdata(1 : length(sdata)/4)));
        Force2p = sgolayfilt(sgolayfilt(Force2,1,49),1,49);


        Save = 'N';

        while Save == 'N'
            
            close all
            fig = figure('units','normalized','outerposition',[0 0 1 1]);
            plot(Force2p)
            hold on
            
            clear('locs', 'pks')

            prompt ='Manuel(1) : choix du seuil et des pics / Automatique(2) : Choix du seuil large / Semi-Automatique(3) : Choix du seuil restreint ? [ 1 / 2 / 3 ]';
            M_A = input(prompt,'s');

            if M_A == '1'
                dcm_obj = datacursormode(fig);
                set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on')

                str='N';
                while str == 'N'
                    prompt ='Plot des points? [Y/N]';
                    str = input(prompt,'s');

                    if str ~= 'N'
                        c_info = getCursorInfo(dcm_obj);
                        hold on
                        pks = Force2p([c_info.DataIndex]);
                        locs = [c_info.DataIndex];
                        plot(locs, pks, 'o')
                        seuil = input('seuil ? ');
                        Force2_pics = Force2p;
                        Force2_pics(Force2_pics<seuil) = 0;

                            for i = 1 : length(pks)
                                dif = Force2_pics(locs(i));
                                while dif > delta
                                    locs(i) = locs(i) - 2 ; 
                                    dif = Force2_pics(locs(i));
                                end
                            end
                            plot(locs , Force2(locs), 'ks')
                    end
                end
            end

            if M_A == '3'
                seuil = input('Valeur du Seuil ? ');
                line(xlim, [seuil seuil], 'Linestyle', '--', 'Color', [1 0 0])

                Force2_pics = Force2p;
                Force2_pics(Force2_pics<seuil) = 0; 

                [pks, locs] = findpeaks(Force2_pics);
                plot(locs, pks, 'o')

                delta = 0.1;

                for i = 1 : length(pks)

                    dif = Force2_pics(locs(i));
                    while dif > delta
                        locs(i) = locs(i) - 2 ; 
                        dif = Force2_pics(locs(i));
                    end
                end
                plot(locs , Force2(locs), 'ks')
            end
            
            if M_A == '2'
                seuil = input('Valeur du Seuil ? ');
                line(xlim, [seuil seuil], 'Linestyle', '--', 'Color', [1 0 0])

                Force2_pics = Force2p;
                Force2_pics(Force2_pics<seuil) = 0; 

                [pks, locs] = findpeaks(Force2_pics);
                plot(locs, pks, 'o')

                delta = 0.1;

                for i = 1 : length(pks)

                    dif = Force2p(locs(i));
                    while dif > delta
                        locs(i) = locs(i) - 2 ; 
                        dif = Force2p(locs(i));
                    end
                end
                plot(locs , Force2(locs), 'ks')
            end
            
            Save = input('Save ? [ Y / N ] ', 's');

            if Save ~= 'N'                
                GRPs(grp).EMG(2).reps2 = locs(find([100 diff(locs)]>99));
            end
        end
    end
%    save([pathname sujet{suj} '\' 'NewREPS.mat'],'GRPs')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Données cinématiques                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for suj = 3

    load([pathname_Cine 'allDataPiano_' sujet{suj}])
       
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    close all 
    
    for grp = 14 %1:16
        
%        fich=mcreadc3d(['\\10.89.24.15\e\Projet_piano\Bassin\' sujet{suj} '\EnregistrementLA\' allREPS(11).sujet(grp).name  '.c3d']);
        
        
        % Force 1, premier serie d'essai 
        a = allData(suj,grp).Tobs(:,68,:);
        data = a(:,:);

%         for j = 1 : length(fich.markerName)
%             index(j) ={fich.markerName{j,1}(end-1:end)};
%         end
        
%         i_touche = find( strcmp(index, 'A4')) ;  

%         data = fich.data( :, i_touche*3-2 : i_touche*3)';
        sdata = sort(data); 
        Cine = -(data(3,:));
        Cine([1:200, length(Cine)-300:length(Cine)]) = NaN;


        Save = 'N';
        while Save == 'N'
            close all
            fig = figure('units','normalized','outerposition',[0 0 1 1]);
            plot(Cine, 'k')
            hold on
            clear('locs', 'pks')
            
            
%            touche = -97.5;
            touche = -0.107; 
            b = find(Cine<touche);
            
            Cine_pics = Cine;
            
            interv = [b(1) b(length(b))];
            Cine_pics([1:interv(1), interv(2):end]) = touche;
            Cine_pics(isnan(Cine_pics))=touche;
            plot( Cine_pics)
            
            line(xlim, [touche touche], 'Color', [0.5 0 1])
            
            decoupe = 'N';
            decoupe = input('decoupe OK ? [Y / N] ', 's');
            
            while decoupe == 'N' 
                
                interv = input('interval d analyse  [a b] ');
                Cine_pics([1:interv(1), interv(2):end]) = touche;
                hold off
                plot(Cine)
                hold on
                plot( Cine_pics)
                line(xlim, [touche touche], 'Color', [0.5 0 1])
                decoupe = input('decoupe OK ? [Y / N] ', 's');
            end
            
            z_touche = 'N';
            z_touche = input('Touche ok ?  [Y / N]', 's');
            
            while z_touche == 'N'
                
                touche = input('valeur touche ? ');
                hold off
                plot(Cine)
                hold on
                plot( Cine_pics)
                clear('locs', 'pks')         
                line(xlim, [touche touche], 'Color', [0.5 0 1])
                
                z_touche = input('Touche ok ?  [Y / N]', 's');
            end
            
            xlim([interv(1), interv(2)])
            prompt = 'Manuel(1) : choix du seuil et des pics / Automatique(2) : Choix du seuil large [ 1 / 2 ]';
            M_A = input(prompt,'s');

            if M_A == '1'
                dcm_obj = datacursormode(fig);
                set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on')

                str='N';
                while str == 'N'
                    prompt ='Plot des points? [Y/N]';
                    str = input(prompt,'s');

                    if str ~= 'N'
                        c_info = getCursorInfo(dcm_obj);
                        hold on
                        pks = Cine([c_info.DataIndex]);
                        locs = [c_info.DataIndex];
                        plot(locs, pks, 'o')
                        
                        seuil = input('seuil ? ');
                        line(xlim, [seuil seuil], 'Linestyle', '--', 'Color', [1 0 0])
                        Cine_pics(Cine_pics<seuil) = touche; 
                        
                            for i = 1 : length(pks)
                                dif = Cine(locs(i));
                                while dif > touche
                                    locs(i) = locs(i) - 2 ; 
                                    dif = Cine(locs(i));
                                end
                            end
                             plot(locs , touche * ones(1 , length(locs)), 'ks')
                    end
                end
            end

            
            if M_A == '2'
                
            seuil = input('Valeur du Seuil ? ');
            line(xlim, [seuil seuil], 'Linestyle', '--', 'Color', [1 0 0])
            
            Cine_pics(Cine_pics<seuil) = touche; 
            
            plot(Cine_pics)

            [pks, locs] = findpeaks(Cine_pics);
            plot(locs, pks, 'o')


            for i = 1 : length(pks)

                dif = Cine(locs(i));
                while dif > touche
                    locs(i) = locs(i) - 2 ; 
                    dif = Cine(locs(i));
                end
            end
            plot(locs , touche * ones(1 , length(locs)), 'ks')
            end

            Save = input('Save ? [ Y / N ] ', 's');

%             if Save ~= 'N'
%                 allREPS(suj).sujet(grp).name = allData(suj,grp).trialName;
%                 allREPS(suj).sujet(grp).rate.cine = 150;
%                 allREPS(suj).sujet(grp).reps_cine = locs(find([100 diff(locs)]>99));
%             end

            if Save ~= 'N'
                allREPS(suj).sujet(grp).reps_cine = locs(find([100 diff(locs)]>99));
            end
        end
    end
    
    
%    load([pathname sujet{suj} '\' 'NewREPS.mat'])
    
%     n = {allREPS(suj).sujet.name};
%     
%     for grp = 1 : 8
%         for i = 1 : 2
%             place = find(strcmp(n,GRPs(grp).EMG(i).Filename));
%             allREPS(suj).sujet(place).reps_fo = GRPs(grp).EMG(i).reps2 ; 
%             allREPS(suj).sujet(place).rate.force = GRPs(grp).EMG(i).Rate ;
%         end
%     end
%     
%     
%     load('\\10.89.24.15\e\Projet_piano\Bassin\Son\Attack_frames_son.mat')
%     
%     for grp = 1 : 16
%         allREPS(suj).sujet(grp).reps_son = Attack_frames_Son(suj).frame(grp, 2:end)*44100;
%         allREPS(suj).sujet(grp).rate.son = 44100 ;
%     end
    
%    save([pathname 'allREPS.mat'],'allREPS')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          Données sonores                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for suj = [1 2 4 5 7 11]
    for grp = 1 : 16
        allREPS(suj).sujet(grp).reps_son = Attack_frames_Son(suj).frame(grp, 2:end)*44100;
        allREPS(suj).sujet(grp).rate.son = 44100 ;
    end
end

%%
for suj = [1 2 4 7]
    for grp = 1:16 

        A = allREPS(suj).sujet(grp).reps_cine   / allREPS(suj).sujet(grp).rate.cine   ;
        B = allREPS(suj).sujet(grp).reps_fo     / allREPS(suj).sujet(grp).rate.force   ;
        C = allREPS(suj).sujet(grp).reps_son    / allREPS(suj).sujet(grp).rate.son   ;
%        D = allREPS(suj).sujet(grp).reps_cine_touche   / allREPS(suj).sujet(grp).rate.cine   ;
        C(C==0)=NaN;
        C = (C(~isnan(C)) - (C(1)- A(1)));    

        %[X,Y,D(grp)] = alignsignals(diff(C(~isnan(C))), diff(A));
        %D = 1;

        [a_l b_l] = max([ length(A) length(B) (length(C)-1)]);

        allREPS(suj).sujet(grp).diff = NaN(3,a_l);

        if b_l == 1

            allREPS(suj).sujet(grp).diff(b_l,:) = A ;

            for i = 1 : length(B)
                [difmin, difloc] = min( abs( diff( [A; B(i)*ones(1,length(A))] )));
                allREPS(suj).sujet(grp).diff(2,difloc) = B(i);
            end

            for i = 1 : length(C)
                [difmin, difloc] = min( abs( diff( [A; C(i)*ones(1,length(A))] )));
                allREPS(suj).sujet(grp).diff(3,difloc) = C(i);
            end
            
%             for i = 1 : length(D)
%                 [difmin, difloc] = min( abs( diff( [A; D(i)*ones(1,length(A))] )));
%                 allREPS(suj).sujet(grp).diff(4,difloc) = D(i);
%             end

        end

        if b_l == 2

            allREPS(suj).sujet(grp).diff(b_l,:) = B ;

            for i = 1 : length(A)
                [difmin, difloc] = min( abs( diff( [B; A(i)*ones(1,length(B))] )));
                allREPS(suj).sujet(grp).diff(1,difloc) = A(i);
            end

            for i = 1 : length(C)
                [difmin, difloc] = min( abs( diff( [B; C(i)*ones(1,length(B))] )));
                allREPS(suj).sujet(grp).diff(3,difloc) = C(i);
            end
            
%             for i = 1 : length(D)
%                 [difmin, difloc] = min( abs( diff( [B; D(i)*ones(1,length(B))] )));
%                 allREPS(suj).sujet(grp).diff(4,difloc) = D(i);
%             end

        end

        if b_l == 3

            allREPS(suj).sujet(grp).diff(b_l,:) = C ;

            for i = 1 : length(A)
                [difmin, difloc] = min( abs( diff( [C; A(i)*ones(1,length(C))] )));
                allREPS(suj).sujet(grp).diff(1,difloc) = A(i);
            end

            for i = 1 : length(B)
                [difmin, difloc] = min( abs( diff( [C; B(i)*ones(1,length(C))] )));
                allREPS(suj).sujet(grp).diff(2,difloc) = B(i);
            end
%             
%             for i = 1 : length(D)
%                 [difmin, difloc] = min( abs( diff( [C; D(i)*ones(1,length(C))] )));
%                 allREPS(suj).sujet(grp).diff(4,difloc) = D(i);
%             end

        end
    end
end

%% 
for suj = [1 2 4 7]
    for grp = 1 : 16
        figure(suj)
        subplot(8,2,grp)
        boxplot(allREPS(1).sujet(grp).diff  )

        m_dif_fo(suj,grp) = nanmean( diff( allREPS(suj).sujet(grp).diff(1:2,:)));
        m_dif_son(suj,grp) = nanmean( diff( allREPS(suj).sujet(grp).diff([1 3],:)));
    end
end

%% Aligner une attaque cinematique avec les données Fo puis caler le son dessus

for suj = 1:12
    
    close all
    
    for grp = 1 : 16
        Cine_rep = allREPS(suj).sujet(grp).reps_cine(1) /allREPS(suj).sujet(grp).rate.cine;
        Fo = allREPS(suj).sujet(grp).reps_fo /allREPS(suj).sujet(grp).rate.force;
        Son = allREPS(suj).sujet(grp).reps_son /allREPS(suj).sujet(grp).rate.son;
        Son(Son == 0) = [];
        
        figure(grp)
        %subplot(4,4,grp)
        plot(Son, 'bo')
        hold on; plot(Fo, 'ro')
        
        [a b] = min(abs(Fo - Cine_rep));
        
        Son = Son - (Son(b) - Cine_rep);
        
%         disB = [0 diff(Fo)];
%         disC = [0 diff(Son)];
%         [difmin, difloc] = min( abs( diff( [disB; disC(b)*ones(1,length(B))] )));
        
%         allREPS(suj).sujet(grp).true_REPS = C - (C(difloc - (b - difloc)) - A);
        allREPS(suj).sujet(grp).true_REPS = Son;
        
       
        hold on
        plot(Cine_rep, 'go')
        plot(Fo, 'ko')
        plot(Son, 'rX', 'Markersize', 12)
        
    end
end
