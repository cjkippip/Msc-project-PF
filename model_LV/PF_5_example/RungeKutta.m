function [dx1, dx2, dx3, dx4] = RungeKutta(x)
% Fourth order Runge Kutta integration for the falling body system.
global rho0 g k dt
dx1(1,1) = x(2);
dx1(2,1) = rho0 * exp(-x(1)/k) * x(2)^2 / 2 * x(3) - g;
dx1(3,1) = 0;
dx1 = dx1 * dt;
xtemp = x + dx1 / 2;
dx2(1,1) = xtemp(2);
dx2(2,1) = rho0 * exp(-xtemp(1)/k) * xtemp(2)^2 / 2 * xtemp(3) - g;
dx2(3,1) = 0;
dx2 = dx2 * dt;
xtemp = x + dx2 / 2;
dx3(1,1) = xtemp(2);
dx3(2,1) = rho0 * exp(-xtemp(1)/k) * xtemp(2)^2 / 2 * xtemp(3) - g;
dx3(3,1) = 0;
dx3 = dx3 * dt;
xtemp = x + dx3;
dx4(1,1) = xtemp(2);
dx4(2,1) = rho0 * exp(-xtemp(1)/k) * xtemp(2)^2 / 2 * xtemp(3) - g;
dx4(3,1) = 0;
dx4 = dx4 * dt;
return;


% function [dx1, dx2, dx3, dx4] = RungeKutta(x)
% % Fourth order Runge Kutta integration for the falling body system.
% global rho0 g k dt
% dx1(1,1) = x(2);
% dx1(2,1) = rho0 * exp(-x(1)/k) * x(2)^2 / 2 * x(3) - g;
% dx1(3,1) = 0;
% dx1 = dx1 * dt;
% xtemp = x + dx1 / 2;
% dx2(1,1) = xtemp(2);
% dx2(2,1) = rho0 * exp(-xtemp(1)/k) * xtemp(2)^2 / 2 * xtemp(3) - g;
% dx2(3,1) = 0;
% dx2 = dx2 * dt;
% xtemp = x + dx2 / 2;
% dx3(1,1) = xtemp(2);
% dx3(2,1) = rho0 * exp(-xtemp(1)/k) * xtemp(2)^2 / 2 * xtemp(3) - g;
% dx3(3,1) = 0;
% dx3 = dx3 * dt;
% xtemp = x + dx3;
% dx4(1,1) = xtemp(2);
% dx4(2,1) = rho0 * exp(-xtemp(1)/k) * xtemp(2)^2 / 2 * xtemp(3) - g;
% dx4(3,1) = 0;
% dx4 = dx4 * dt;
% return;








