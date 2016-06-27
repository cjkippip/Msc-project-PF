figure(1),clf,
plot(xx1,BSOptionCPrice,'b','LineWidth',1.5);
axis([56,222,0,inf])
title('Comparison between BS and true price of call option','FontSize',15)
ylabel('Call option price','FontSize',13,'FontWeight','bold')
xlabel('time(T/4+1 to T)','FontSize',13,'FontWeight','bold')
hold on
plot(xx1,trueOptionCPrice,'r','LineWidth',1.5);
legend({'BS price','True price'},'Location','northeast','FontSize',13,'FontWeight','bold');
grid on
grid minor
hold off
%%
a=wgn(2, 1, 10*log10(5));
%%
nIx1 = dims~=2;
nIx2 = dims~=1;
aa=[1,2];
aa1=aa(nIx1);
aa2=aa(nIx2);
%%
figure(7),clf,
scatter(PCenter(1,1),PCenter(2,1),'r.','LineWidth ',1.5);


