%% plot tureState
figure(10),clf,
plot(trueState(1,:),'r')
hold on
plot(trueState(2,:),'b')
plot(trueState(3,:),'g')
grid on
%%
a1=2;
a2=[1 2 3;4 5 6;7 8 9];
a3=a1*a2;
%%
a1=randn([2,10000000]);
a2=mean(a1,2);



