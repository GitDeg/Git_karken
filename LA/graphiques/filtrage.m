F=2100;

t = 1/F-525/F: 1/F : 525/F;

figure
subplot (4,1,1)
plot(t , GRP(2).EMG(4).cycle(11,:)')
legend('unfiltred cycle')
title('A')

subplot(4,1,2)
plot(t , GRP(2).EMG(4).cycle_filt(11,:)')
legend('bandpass filter')
title('B')

subplot(4,1,3)

plot(t , GRP(2).EMGal(4).cycle_ali_env(11,:)')
hold on
plot(t , abs(GRP(2).EMG(4).cycle_filt(11,:)'),'color',[0,0,1,0.3])
legend('RMS sliding window','rectified cycle')

title('C')

subplot(4,1,4)
hold on
plot(t , GRP(2).EMGal(4).cycle_ali_env(1,:),'color','k')
plot(t , GRP(2).EMGal(4).cycle_ali_env','color','k','HandleVisibility','off')

plot(t(1:end-1) ,GRP(2).EMGal(4).cycle_ali_ttt(1,:),'color','b')
plot(t(1:end-1) ,GRP(2).EMGal(4).cycle_ali_ttt','color','b','HandleVisibility','off')
plot(t(1:end-1) ,GRP(2).EMGal(4).mean_ttt','color','r','linewidth',2)
legend('All cycles','Statistical filtering','Mean EMG')

title('D')
xlabel('Temps (s)')