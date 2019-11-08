figure
subplot(4,1,1)
plot(All_data(1).CINE(1).mark(5).data((1:10)*3-2,:)')
ylim([-75 150])
subplot(4,1,2)
plot(All_data(2).CINE(1).mark(5).data((1:10)*3-2,:)')
ylim([-75 150])
subplot(4,1,3)
plot(All_data(5).CINE(1).mark(1).data((1:10)*3-2,:)')
ylim([-75 150])
subplot(4,1,4)
plot(All_data(6).CINE(1).mark(1).data((1:10)*3-2,:)')
ylim([-75 150])

%%
for grp = 1 : 6
    for suj = 1 : 12
        for mu = 1 : 9
            
            d(mu,suj,grp) = 1 / (mean(max(All_data(grp).EMG(suj).mu(mu).data_MVC'))- mean(min(All_data(grp).EMG(suj).mu(mu).data_MVC')));
            
        end
    end
end
figure
subplot(3,2,1)
boxplot(d(ordre_mu,:,1)')
subplot(3,2,2)
boxplot(d(ordre_mu,:,2)')
subplot(3,2,3)
boxplot(d(ordre_mu,:,3)')
subplot(3,2,4)
boxplot(d(ordre_mu,:,4)')
subplot(3,2,5)
boxplot(d(ordre_mu,:,5)')
subplot(3,2,6)
boxplot(d(ordre_mu,:,6)')

%%

figure
subplot(4,1,1)
plot(All_data(1).EMG(1).mu(2).data_MVC')
ylim([0 0.1])
subplot(4,1,2)
plot(All_data(2).EMG(1).mu(2).data_MVC')
ylim([0 0.1])
subplot(4,1,3)
plot(All_data(5).EMG(1).mu(2).data_MVC')
ylim([0 0.1])
subplot(4,1,4)
plot(All_data(6).EMG(1).mu(2).data_MVC')
ylim([0 0.1])
%%
figure
subplot(2,1,1)
plot(All_data(1).CINE(1).mark(5).data'  )
subplot(2,1,2)
plot(All_data(5).CINE(1).mark(1).data'  )



figure
A = All_data(1).CINE(1).mark(5).data((1:10)*3,:);
B = All_data(5).CINE(1).mark(1).data((1:10)*3,:);
TFA = islocalmin(A,2);
TFB = islocalmin(B,2);

hold on 
for i = 1 : 10
    plot( All_data(1).CINE(1).mark(5).data(i*3-2,TFA(i,:)),  All_data(1).CINE(1).mark(5).data(i*3-1,TFA(i,:)), 'bX')
    plot( All_data(5).CINE(1).mark(1).data(i*3-2,TFB(i,:)),  All_data(5).CINE(1).mark(1).data(i*3-1,TFB(i,:)), 'rO')
end
%%
for grp = 1 : 6
    for suj = 1 : 12
        if grp < 5 
            mark = 5;
        else
            mark = 1;
        end
        TF = islocalmin(All_data(grp).CINE(suj).mark(mark).data((1:10)*3,:), 2);
        
        All_data(grp).CINE(suj).end_p = [];
        for i = 1 : 10
            if grp < 5 && sum(TF(i,:) == 1) < 4
                TF(i,end) = 1;
            end
            All_data(grp).CINE(suj).end_p = [All_data(grp).CINE(suj).end_p , [All_data(grp).CINE(suj).mark(mark).data(i*3-2,TF(i,:));...
                                                                                All_data(grp).CINE(suj).mark(mark).data(i*3-1,TF(i,:))]];
        end
    end
end

