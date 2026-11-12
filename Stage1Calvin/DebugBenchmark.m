
clear; clc; close all;

%% Parameters
theta = [1.5, 0.5, 1.0, 0.5];

t0 = 0;
tEnd = 100;        % Integration endpoint
u0 = [2.11; 1.35];

hMax = 0.01;       % RK4 step size

%% Lotka-Volterra system
f = @(t,u) [
    theta(1)*u(1) - theta(2)*u(1)*u(2);
   -theta(3)*u(2) + theta(4)*u(1)*u(2)
];

%% ==========================
%  Classical RK4
%  ==========================

N = ceil((tEnd - t0)/hMax);
h = (tEnd - t0)/N;

tRK4 = linspace(t0, tEnd, N+1)';
URK4 = zeros(N+1, 2);
URK4(1,:) = u0';

tic;

for n = 1:N

    u = URK4(n,:)';
    t = tRK4(n);

    k1 = f(t, u);
    k2 = f(t+h/2, u+h*k1/2);
    k3 = f(t+h/2, u+h*k2/2);
    k4 = f(t+h, u+h*k3);

    URK4(n+1,:) = ...
        (u + h*(k1+2*k2+2*k3+k4)/6)';

end

timeRK4 = toc;

%% ==========================
%  ODE45
%  ==========================

opts = odeset( ...
    'RelTol', 1e-11, ...
    'AbsTol', 1e-13);

tic;

% Two endpoints: retain adaptive solver history
sol = ode45(f, [t0, tEnd], u0, opts);

timeODE45 = toc;

% Every accepted ODE45 step
tODE45 = sol.x(:);
UODE45 = sol.y';

%% ==========================
%  Compare at identical times
%  ==========================

% Evaluate ODE45's continuous solution at RK4 times
UODE45_interp = deval(sol, tRK4')';

errorX = abs(URK4(:,1) - UODE45_interp(:,1));
errorY = abs(URK4(:,2) - UODE45_interp(:,2));

fprintf('\n===== Solver Comparison =====\n');
fprintf('RK4 time:   %.6f seconds\n', timeRK4);
fprintf('ODE45 time: %.6f seconds\n', timeODE45);
fprintf('RK4 steps:  %d\n', N);
fprintf('ODE45 accepted steps: %d\n', numel(tODE45)-1);
fprintf('Max error x: %.12e\n', max(errorX));
fprintf('Max error y: %.12e\n', max(errorY));

%% ==========================
%  Plot complete solver history
%  ==========================

figure('Color','w');

subplot(3,1,1);
plot(tRK4, URK4(:,1), 'b-', 'LineWidth', 1);
hold on;
plot(tODE45, UODE45(:,1), 'r--', 'LineWidth', 1);
xlabel('Time');
ylabel('Prey x');
legend('RK4','ODE45');
title('Complete Prey Trajectory');
grid on;

subplot(3,1,2);
plot(tRK4, URK4(:,2), 'b-', 'LineWidth', 1);
hold on;
plot(tODE45, UODE45(:,2), 'r--', 'LineWidth', 1);
xlabel('Time');
ylabel('Predator y');
legend('RK4','ODE45');
title('Complete Predator Trajectory');
grid on;

subplot(3,1,3);
semilogy(tRK4, max(errorX,eps), 'b-');
hold on;
semilogy(tRK4, max(errorY,eps), 'r-');
xlabel('Time');
ylabel('Absolute Difference');
legend('Error x','Error y');
title('RK4 vs ODE45');
grid on;

%% ==========================
%  Phase portrait
%  ==========================

figure('Color','w');

plot(URK4(:,1), URK4(:,2), 'b-', 'LineWidth', 1);
hold on;
plot(UODE45(:,1), UODE45(:,2), 'r--', 'LineWidth', 1);
plot(u0(1), u0(2), 'ko', 'MarkerFaceColor','k');

xlabel('Prey x');
ylabel('Predator y');
legend('RK4','ODE45','Initial Condition');
title('Lotka-Volterra Phase Portrait');
grid on;

%% ==========================
%  Conservation law diagnostic
%  ==========================

% H(x,y) = theta4*x - theta3*log(x)
%        + theta2*y - theta1*log(y)

H = @(U) theta(4)*U(:,1) ...
       - theta(3)*log(U(:,1)) ...
       + theta(2)*U(:,2) ...
       - theta(1)*log(U(:,2));

HRK4 = H(URK4);
HODE = H(UODE45);

figure('Color','w');

plot(tRK4, HRK4-HRK4(1), 'b-', 'LineWidth', 1);
hold on;
plot(tODE45, HODE-HODE(1), 'r--', 'LineWidth', 1);

xlabel('Time');
ylabel('H(t) - H(0)');
legend('RK4','ODE45');
title('Conservation Law Error');
grid on;

fprintf('\n===== Conservation Law =====\n');
fprintf('RK4 max drift:   %.12e\n', ...
    max(abs(HRK4-HRK4(1))));
fprintf('ODE45 max drift: %.12e\n', ...
    max(abs(HODE-HODE(1))));
