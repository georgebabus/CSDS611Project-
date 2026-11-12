
function [x, y] = LV_ODE45(theta, tObs, u0)

% Lotka-Volterra right-hand side
f = @(t, u) [
    theta(1)*u(1) - theta(2)*u(1)*u(2);
    -theta(3)*u(2) + theta(4)*u(1)*u(2)
    ];

% Solver options
opts = odeset('RelTol', 1e-10, ...
    'AbsTol', 1e-12);

% Solve at observation times
[~, U] = ode45(f, tObs, u0(:), opts);

% Extract prey and predator populations
x = U(:, 1);
y = U(:, 2);

end
