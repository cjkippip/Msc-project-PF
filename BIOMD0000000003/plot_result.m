%% Plot result

%% estimation of states and parameter

% state
% figure(3),clf,
% for k=1:state_dim                               
%     subplot(state_dim,1,k)
%     plot(tspan, trueState(k,:), 'b-', tspan, meanP_State(k,:), 'r-')
%     if k==1
%         title('PF estimate state','FontSize',15);
%     end
%     xlabel('Time','FontSize',13,'FontWeight','bold'); 
%     ylabel('Number of spicies','FontSize',13,'FontWeight','bold'); 
%     legend({'ture','estimation'},...
%         'Location','eastoutside','FontSize',11,'FontWeight','bold'); 
%     grid on
%     grid minor
% end

% parrameter
figure(4),clf,
for k=1:3    
    subplot(3,1,k)
    
    plot(tspan, truePara(k,:), 'b--', tspan, meanP_Para(k,:), 'r-')
    
    xlabel('Time','FontSize',13,'FontWeight','bold');    
    ylabel(['Parameter',allParaName(k)],...
    'FontSize',13,'FontWeight','bold'); 
    
    legend({'ture','estimation'},'Location','eastoutside',...
        'FontSize',11,'FontWeight','bold'); 
    grid on
    grid minor
end
figure(5),clf,
for k=1:3
    index_para_show=k+3;
    subplot(3,1,k)
    
    plot(tspan, truePara(index_para_show,:),...
        'b--', tspan, meanP_Para(index_para_show,:), 'r-')
    
    xlabel('Time','FontSize',13,'FontWeight','bold'); 
    ylabel(['Parameter',allParaName(index_para_show)],...
        'FontSize',13,'FontWeight','bold');
    
    legend({'ture','estimation'},'Location','eastoutside',...
        'FontSize',11,'FontWeight','bold');  
    grid on
    grid minor
end

%% plot distribution of states and parameters
% state
% figure(6),clf,
% subplot(211)
% hist(P_StateSpace(1,:,T),50);
% hold on
% plot(meanP_State(1,T),0,'r*');
% 
% subplot(212)
% hist(P_StateSpace(2,:,T),50);
% hold on
% plot(meanP_State(2,T),0,'r*');

% parameter
figure(7),clf,
subplot(231)
hist(P_ParaSpace(1,:,T),50);
hold on
plot(meanP_Para(1,T),0,'r*');
plot(alphaStar,0,'ro');
title('true k6 = 1.0')

subplot(232)
hist(P_ParaSpace(2,:,T),50);
hold on
plot(meanP_Para(2,T),0,'r*');
plot(betaStar,0,'ro');
title('true k8notP = 100000.0')

subplot(233)
hist(P_ParaSpace(3,:,T),50);
hold on
plot(meanP_Para(3,T),0,'r*');
plot(deltaStar,0,'ro');
title('true k9 = 1000.0')

subplot(234)
hist(P_ParaSpace(4,:,T),50);
hold on
plot(meanP_Para(4,T),0,'r*');
plot(gammaStar,0,'ro');
title('true k3 = 200.0')

subplot(235)
hist(P_ParaSpace(5,:,T),50);
hold on
plot(meanP_Para(5,T),0,'r*');
plot(gammaStar,0,'ro');
title('true k5notP = 0.0')

subplot(236)
hist(P_ParaSpace(6,:,T),50);
hold on
plot(meanP_Para(6,T),0,'r*');
plot(gammaStar,0,'ro');
title('true k1aa = 0.015')

figure(8),clf,
subplot(221)
hist(P_ParaSpace(7,:,T),50);
hold on
plot(meanP_Para(7,T),0,'r*');
plot(gammaStar,0,'ro');
title('true k2 = 0.0')

subplot(222)
hist(P_ParaSpace(8,:,T),50);
hold on
plot(meanP_Para(8,T),0,'r*');
plot(gammaStar,0,'ro');
title('true k7 = 0.6')

subplot(223)
hist(P_ParaSpace(9,:,T),50);
hold on
plot(meanP_Para(9,T),0,'r*');
plot(gammaStar,0,'ro');
title('true k4 = 180.0')

subplot(224)
hist(P_ParaSpace(10,:,T),50);
hold on
plot(meanP_Para(10,T),0,'r*');
plot(gammaStar,0,'ro');
title('true k4prime = 0.018')


