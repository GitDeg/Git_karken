load('allREPS.mat')

GrpNames = {' BasFraSta ', ...
            ' BasFraTen ', ...
            ' BasPreSta ', ...
            ' BasPreTen ', ...
            ' MsFraSta ', ...
            ' MsFraTen ', ...
            ' MsPreSta ', ...
            ' MsPreTen '};
            
for grp = 1 : 16
    
    close all
    i = 0;
    for suj = [2 4 7 11]
        i = i+1;
        
        Cine = allREPS(suj).sujet(grp).diff(1,1:20)  ;
        Fo = allREPS(suj).sujet(grp).diff(2,1:20)  ;
        
        plot(Cine, 'rx')
        hold on ;
        plot(Fo, 'bo')
        
        dif_cine_fo(i,:) = diff([Cine; Fo]); % Fo - Cine donc si > 0, Cine avant Force.  Si < 0, Force avant Cine. 
    end
    
    boxplot_data(:,grp) = dif_cine_fo(:) ;
    boxplot(boxplot_data, [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8], 'labels', GrpNames)

end

y = boxplot_data(:);

g1 = [repmat({'Bas'},160,1); repmat({'Bas'},160,1); repmat({'Bas'},160,1); repmat({'Bas'},160,1); repmat({'Ms'},160,1); repmat({'Ms'},160,1); repmat({'Ms'},160,1); repmat({'Ms'},160,1)] ;
g2 = [repmat({'Fra'},160,1); repmat({'Fra'},160,1); repmat({'Pre'},160,1); repmat({'Pre'},160,1); repmat({'Fra'},160,1); repmat({'Fra'},160,1); repmat({'Pre'},160,1); repmat({'Pre'},160,1)] ;
g3 = [repmat({'Sta'},160,1); repmat({'Ten'},160,1); repmat({'Sta'},160,1); repmat({'Ten'},160,1); repmat({'Sta'},160,1); repmat({'Ten'},160,1); repmat({'Sta'},160,1); repmat({'Ten'},160,1)] ;

p = anovan(y,{g1,g2,g3});
