figure
plot3(Data(1).CINE(1).data(1,:), Data(1).CINE(1).data(2,:), Data(1).CINE(1).data(3,:))
hold on 
plot3(Data(1).CINE(2).data(1,:), Data(1).CINE(2).data(2,:), Data(1).CINE(2).data(3,:))

plot3(CINE_SAUT(11).RawData(41).Data(1,:), CINE_SAUT(11).RawData(41).Data(2,:), CINE_SAUT(11).RawData(41).Data(3,:))

plot3(CINE_SAUT(19).RawData(41).Data(1,:), CINE_SAUT(19).RawData(41).Data(2,:), CINE_SAUT(19).RawData(41).Data(3,:))
grid on


%%
figure
subplot(3,1,1)
plot(CINE_SAUT(4).RawData(40).Data')

hold on 
subplot(3,1,2)
plot(CINE_SAUT(12).RawData(41).Data')
subplot(3,1,3)
plot(CINE_SAUT(20).RawData(41).Data')
%%
A = CINE_SAUT(4).RawData(40).Data(1,:);
plot(A)

[pk, lc] = findpeaks(A,'MinPeakDistance', 80);

%%
B = EMG_SAUT(4).RawData(9,:);
plot((1:length(EMG_SAUT(4).RawData(9,:)))/2100, B)
yyaxis right
plot((1:length(A))/150, A)

line([reps(1:end-1); reps(1:end-1)], repmat(ylim, length(reps(1:end-1)),1)' )

%%
for i = 2 : length(reps(1:end-1))
    C(i-1).data = B(reps(i-1)*2100:reps(i)*2100);
end

for i = 1 : 19
    hold on
    plot(C(i).data)
end

%%
figure
A = CINE_SAUT(12).RawData(40).Data(1,:);
plot(A)

[pk, lc] = findpeaks(A,'MinPeakDistance', 80);

lc(lc<640)=[]
lc(lc>4150)=[]
reps = lc(1:2:end)/150

%%
figure
B = EMG_SAUT(12).RawData(9,:);
plot((1:length(B))/2100, B)
yyaxis right
plot((1:length(A))/150, A)

line([reps; reps], repmat(ylim, length(reps),1)' )

%%
for i = 2 : length(reps)
    C(i-1).data = B(round(reps(i-1)*2100): round(reps(i)*2100));
end

figure
for i = 1 : 16
    hold on
    plot(C(i).data)
end