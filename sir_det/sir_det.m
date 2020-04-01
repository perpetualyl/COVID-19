% Script to test time-integration of the SIR model

clear
close all

alpha = 2.5;

% Initial condition
y0 = [0.45; 0.35];
y0 = [y0; 1 - sum(y0)];

% Right-hand side
yprime = @(tt, yy) sir_rhs(yy, alpha);

t0 = 0;
dt = 0.01;
N = 500;

% Time integration and output
Y = rk4(y0, yprime, t0, dt, N);

t = 0:dt:(N*dt);
plot(t, Y(1,:), 'r'); hold on;
plot(t, Y(2,:), 'b'); 
plot(t, Y(3,:), 'k'); 
legend('S', 'I', 'R');