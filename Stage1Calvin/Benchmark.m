
clear; clc; close all;

%% Prescribed parameters
theta = [1.5, 0.5, 1.0, 0.5];

%% Load data and prepend initial condition
Data = readtable('LV_data.csv');

InitialRow = table(0, 2.11, 1.35, ...
    'VariableNames', {'t', 'x', 'y'});

Data = [InitialRow; Data];

tObs = Data.t;
u0 = [Data.x(1); Data.y(1)];

%% Settings
hMax = 0.03;          % RK4 maximum step size
runFullData = true;  % Set true for full-data comparison

%% Part 1: Short-time sanity check
tCheck = (0:0.1:500)';

% RK4
tic;
[xRK4, yRK4] = LV_RK4(theta, tCheck, u0, hMax);
timeRK4 = toc;

% ODE45
tic;
[xODE, yODE] = LV_ODE45(theta, tCheck, u0);
timeODE = toc;

% Compare numerical solutions
errorX = max(abs(xRK4 - xODE));
errorY = max(abs(yRK4 - yODE));

fprintf('\n===== Short-Time Sanity Check =====\n');
fprintf('RK4 time:    %.6f seconds\n', timeRK4);
fprintf('ODE45 time:  %.6f seconds\n', timeODE);
fprintf('Max |xRK4 - xODE45|: %.6e\n', errorX);
fprintf('Max |yRK4 - yODE45|: %.6e\n', errorY);

%% Plot short-time comparison
figure;

subplot(2,1,1);
plot(tCheck, xRK4, 'b-', 'LineWidth', 1.5);
hold on;
plot(tCheck, xODE, 'r--', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Prey x');
legend('RK4', 'ODE45');
title('Prey Population Comparison');
grid on;

subplot(2,1,2);
plot(tCheck, yRK4, 'b-', 'LineWidth', 1.5);
hold on;
plot(tCheck, yODE, 'r--', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Predator y');
legend('RK4', 'ODE45');
title('Predator Population Comparison');
grid on;

%% Part 2: Full-data objective comparison
if runFullData

    % RK4
    tic;
    [xRK4, yRK4] = LV_RK4(theta, tObs, u0, hMax);
    timeRK4 = toc;

    % ODE45
    tic;
    [xODE, yODE] = LV_ODE45(theta, tObs, u0);
    timeODE = toc;

    % Least-squares objectives
    ObjectiveRK4 = sum( ...
        (xRK4(2:end) - Data.x(2:end)).^2 + ...
        (yRK4(2:end) - Data.y(2:end)).^2);

    ObjectiveODE = sum( ...
        (xODE(2:end) - Data.x(2:end)).^2 + ...
        (yODE(2:end) - Data.y(2:end)).^2);

    % Compare results
    fprintf('\n===== Full-Data Comparison =====\n');
    fprintf('RK4 objective:   %.12e\n', ObjectiveRK4);
    fprintf('ODE45 objective: %.12e\n', ObjectiveODE);
    fprintf('RK4 time:        %.6f seconds\n', timeRK4);
    fprintf('ODE45 time:      %.6f seconds\n', timeODE);

    fprintf('Max difference x: %.6e\n', ...
        max(abs(xRK4 - xODE)));
    fprintf('Max difference y: %.6e\n', ...
        max(abs(yRK4 - yODE)));

    %% Plot fitted values against measurements
    figure;

    subplot(2,1,1);
    plot(tObs, xRK4, 'b-o');
    hold on;
    plot(tObs, xODE, 'r--x');
    plot(Data.t(2:end), Data.x(2:end), 'k.');
    xlabel('Time');
    ylabel('Prey x');
    legend('RK4', 'ODE45', 'Data');
    grid on;

    subplot(2,1,2);
    plot(tObs, yRK4, 'b-o');
    hold on;
    plot(tObs, yODE, 'r--x');
    plot(Data.t(2:end), Data.y(2:end), 'k.');
    xlabel('Time');
    ylabel('Predator y');
    legend('RK4', 'ODE45', 'Data');
    grid on;

end
