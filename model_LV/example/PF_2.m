N = 1000;   
Q = 5;      
R = 10;      
T = 10;     
theta = pi/T;       
distance = 80/T;    
WorldSize = 100;    
X = zeros(2, T);    
Z = zeros(2, T);    
P = zeros(2, N);    
PCenter = zeros(2, T);  
w = zeros(N, 1);         
err = zeros(1,T);    
X(:, 1) = [50; 20];     
Z(:, 1) = [50; 20] + wgn(2, 1, 10*log10(R));    
for i = 1 : N
    P(:, i) = [WorldSize*rand; WorldSize*rand];
    dist = norm(P(:, i)-Z(:, 1));     
    w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R); 
end
PCenter(:, 1) = sum(P, 2) / N; 
%%
err(1) = norm(X(:, 1) - PCenter(:, 1));
figure(1),clf,
set(gca,'FontSize',12);
hold on
plot(X(1, 1), X(2, 1), 'b.', 'markersize',30)  
axis([0 100 0 100]);
plot(P(1, :), P(2, :), 'k.', 'markersize',5);   
plot(PCenter(1, 1), PCenter(2, 1), 'r.', 'markersize',25); 
legend('True State', 'Particles', 'The Center of Particles');
title('Initial State');
hold off
grid on
grid minor
%%
for k = 2 : T
    X(:, k) = X(:, k-1) + distance * [(-cos(k * theta)); sin(k * theta)] + wgn(2, 1, 10*log10(Q));     
    Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));    
    for i = 1 : N
        P(:, i) = P(:, i) + distance * [-cos(k * theta); sin(k * theta)] + wgn(2, 1, 10*log10(Q));
        dist = norm(P(:, i)-Z(:, k));     
        w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R); 
    end
    wsum = sum(w);
    for i = 1 : N
        w(i) = w(i) / wsum;
    end 
    
    for i = 1 : N
        wmax = 2 * max(w) * rand; 
        index = randi(N, 1);
        while(wmax > w(index))
            wmax = wmax - w(index);
            index = index + 1;
            if index > N
                index = 1;
            end          
        end
        P(:, i) = P(:, index);     
    end
    
    PCenter(:, k) = sum(P, 2) / N;    
    err(k) = norm(X(:, k) - PCenter(:, k));   
    figure(2),clf,
    set(gca,'FontSize',12);
    clf;
    hold on
    plot(X(1, k), X(2, k), 'b.', 'markersize',50); 
    axis([0 100 0 100]);
    plot(P(1, :), P(2, :), 'k.', 'markersize',5);  
    plot(PCenter(1, k), PCenter(2, k), 'r.', 'markersize',25);
    legend('True State', 'Particle', 'The Center of Particles','Location','southwest');
    hold off
    pause(0.1);
end
 
%%
figure(3),clf,
set(gca,'FontSize',12);
plot(X(1,:), X(2,:), 'b', Z(1,:), Z(2,:), 'g', PCenter(1,:), PCenter(2,:), 'r-','LineWidth',1.5);
hold on
plot(PCenter(1,1),PCenter(2,1),'r.','markersize',25);
axis([0 100 0 100]);
legend({'True State', 'Measurement', 'Particle Filter'},...
    'Location','northwest','FontSize',11,'FontWeight','bold');
xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20);
title('Arc motion', 'FontSize', 15);
hold off
grid on
grid minor
%%
figure(4),clf
set(gca,'FontSize',12);
plot(err,'.-');
xlabel('t=1:10', 'FontSize', 13);
ylabel('Error', 'FontSize', 13);
title('Error(Euclidean distance)', 'FontSize', 15);
grid on
grid minor

