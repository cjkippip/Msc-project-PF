% function [StdErr, EKFErr] = ParticleEx5

% EKF Particle filter example.
% Track a body falling through the atmosphere.
% This system is taken from [Jul00], which was based on [Ath68].
% Compare the particle filter with the EKF particle filter.

global rho0 g k dt

rho0 = 2; % lb-sec^2/ft^4
g = 32.2; % ft/sec^2
k = 2e4; % ft
R = 10^4; % measurement noise variance (ft^2)
Q = diag([0 0 0]); % process noise covariance
M = 10^5; % horizontal range of position sensor
a = 10^5; % altitude of position sensor
P = diag([1e6 4e6 10]); % initial estimation error covariance

x = [3e5; -2e4; 1e-3]; % initial state
xhat = [3e5; -2e4; 1e-3]; % initial state estimate

N = 200; % number of particles  

% Initialize the particle filter.
for i = 1 : N
    xhatplus(:,i) = x + sqrt(P) * [randn; randn; randn]; % standard particle filter
    xhatplusEKF(:,i) = xhatplus(:,i); % EKF particle filter
    Pplus(:,:,i) = P; % initial EKF particle filter estimation error covariance
end

T = 0.5; % measurement time step
randn('state',sum(100*clock)); % random number generator seed

tf = 30; % simulation length (seconds)
dt = 0.5; % time step for integration (seconds)
xArray = x;
xhatArray = xhat;
xhatEKFArray = xhat;

for t = T : T : tf
    fprintf('.');
    % Simulate the system.
    for tau = dt : dt : T
        % Fourth order Runge Kutta ingegration
        [dx1, dx2, dx3, dx4] = RungeKutta(x);
        x = x + (dx1 + 2 * dx2 + 2 * dx3 + dx4) / 6;
        x = x + sqrt(dt * Q) * [randn; randn; randn] * dt;
    end
    % Simulate the noisy measurement.
    z = sqrt(M^2 + (x(1)-a)^2) + sqrt(R) * randn;
    % Simulate the continuous-time part of the particle filters (time update).
    xhatminus = xhatplus;
    xhatminusEKF = xhatplusEKF;
    for i = 1 : N
        for tau = dt : dt : T
            % Fourth order Runge Kutta ingegration
            % standard particle filter
            [dx1, dx2, dx3, dx4] = RungeKutta(xhatminus(:,i));
            xhatminus(:,i) = xhatminus(:,i) + (dx1 + 2 * dx2 + 2 * dx3 + dx4) / 6;
            xhatminus(:,i) = xhatminus(:,i) + sqrt(dt * Q) * [randn; randn; randn] * dt;
            xhatminus(3,i) = max(0, xhatminus(3,i)); % the ballistic coefficient cannot be negative
            % EKF particle filter
            [dx1, dx2, dx3, dx4] = RungeKutta(xhatminusEKF(:,i));
            xhatminusEKF(:,i) = xhatminusEKF(:,i) + (dx1 + 2 * dx2 + 2 * dx3 + dx4) / 6;
            xhatminusEKF(:,i) = xhatminusEKF(:,i) + sqrt(dt * Q) * [randn; randn; randn] * dt;
            xhatminusEKF(3,i) = max(0, xhatminusEKF(3,i)); % the ballistic coefficient cannot be negative
        end
        % standard particle filter
        zhat = sqrt(M^2 + (xhatminus(1,i)-a)^2);
        vhat(i) = z - zhat;
        % EKF particle filter
        zhatEKF = sqrt(M^2 + (xhatminusEKF(1,i)-a)^2);
        F = [0 1 0; -rho0 * exp(-xhatminusEKF(1,i)/k) * xhatminusEKF(2,i)^2 / 2 / k * xhatminusEKF(3,i) ...
            rho0 * exp(-xhatminusEKF(1,i)/k) * xhatminusEKF(2,i) * xhatminusEKF(3,i) ...
            rho0 * exp(-xhatminusEKF(1,i)/k) * xhatminusEKF(2,i)^2 / 2; ...
            0 0 0];
        H = [(xhatminusEKF(1,i) - a) / sqrt(M^2 + (xhatminusEKF(1,i)-a)^2) 0 0];
        Pminus(:,:,i) = F * Pplus(:,:,i) * F' + Q;
        K = Pminus(:,:,i) * H' * inv(H * Pminus(:,:,i) * H' + R);
        xhatminusEKF(:,i) = xhatminusEKF(:,i) + K * (z - zhatEKF);
        zhatEKF = sqrt(M^2 + (xhatminusEKF(1,i)-a)^2);
        vhatEKF(i) = z - zhatEKF;
    end
    % Note that we need to scale all of the q(i) probabilities in a way
    % that does not change their relative magnitudes.
    % Otherwise all of the q(i) elements will be zero because of the
    % large value of the exponential.
    % standard particle filter
    vhatscale = max(abs(vhat)) / 4;
    qsum = 0;
    for i = 1 : N
        q(i) = exp(-(vhat(i)/vhatscale)^2);
        qsum = qsum + q(i);
    end
    % Normalize the likelihood of each a priori estimate.
    for i = 1 : N
        q(i) = q(i) / qsum;
    end
    % EKF particle filter
    vhatscaleEKF = max(abs(vhatEKF)) / 4;
    qsumEKF = 0;
    for i = 1 : N
        qEKF(i) = exp(-(vhatEKF(i)/vhatscaleEKF)^2);
        qsumEKF = qsumEKF + qEKF(i);
    end
    % Normalize the likelihood of each a priori estimate.
    for i = 1 : N
        qEKF(i) = qEKF(i) / qsumEKF;
    end
    % Resample the standard particle filter
    for i = 1 : N
        u = rand; % uniform random number between 0 and 1
        qtempsum = 0;
        for j = 1 : N
            qtempsum = qtempsum + q(j);
            if qtempsum >= u
                xhatplus(:,i) = xhatminus(:,j);
                xhatplus(3,i) = max(0,xhatplus(3,i)); % the ballistic coefficient cannot be negative
                break;
            end
        end
    end
    % The standard particle filter estimate is the mean of the particles.
    xhat = mean(xhatplus')';
    % Resample the EKF particle filter
    Ptemp = Pplus;
    for i = 1 : N
        u = rand; % uniform random number between 0 and 1
        qtempsum = 0;
        for j = 1 : N
            qtempsum = qtempsum + qEKF(j);
            if qtempsum >= u
                xhatplusEKF(:,i) = xhatminusEKF(:,j);
                xhatplusEKF(3,i) = max(0,xhatplusEKF(3,i)); % the ballistic coefficient cannot be negative
                Pplus(:,:,i) = Ptemp(:,:,j);
                break;
            end
        end
    end
    % The EKF particle filter estimate is the mean of the particles.
    xhatEKF = mean(xhatplusEKF')';
    % Save data for plotting.
    xArray = [xArray x];
    xhatArray = [xhatArray xhat];
    xhatEKFArray = [xhatEKFArray xhatEKF];
end

close all;
t = 0 : T : tf;
figure; 
semilogy(t, abs(xArray(1,:) - xhatArray(1,:)), 'b-'); hold;
semilogy(t, abs(xArray(1,:) - xhatEKFArray(1,:)), 'r:'); 
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('Altitude Estimation Error');
legend('Standard particle filter', 'EKF particle filter');

figure; 
semilogy(t, abs(xArray(2,:) - xhatArray(2,:)), 'b-'); hold;
semilogy(t, abs(xArray(2,:) - xhatEKFArray(2,:)), 'r:');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('Velocity Estimation Error');
legend('Standard particle filter', 'EKF particle filter');

figure; 
semilogy(t, abs(xArray(3,:) - xhatArray(3,:)), 'b-'); hold;
semilogy(t, abs(xArray(3,:) - xhatEKFArray(3,:)), 'r:'); 
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('Ballistic Coefficient Estimation Error');
legend('Standard particle filter', 'EKF particle filter');

figure;
plot(t, xArray(1,:));
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('True Position');

figure;
plot(t, xArray(2,:));
title('Falling Body Simulation', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('True Velocity');

for i = 1 : 3
    StdErr(i) = sqrt((norm(xArray(i,:) - xhatArray(i,:)))^2 / tf / dt);
    EKFErr(i) = sqrt((norm(xArray(i,:) - xhatEKFArray(i,:)))^2 / tf / dt);
end
disp(['Standard particle filter RMS error = ', num2str(StdErr)]);
disp(['EKF particle filter RMS error = ', num2str(EKFErr)]);

