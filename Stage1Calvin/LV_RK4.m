
function [x, y] = LV_RK4(theta, tObs, u0, hMax)

    % Initial condition provided as argument
    u = u0(:);

    % Lotka-Volterra right-hand side
    f = @(u) [
        theta(1)*u(1) - theta(2)*u(1)*u(2);
       -theta(3)*u(2) + theta(4)*u(1)*u(2)
    ];

    N = numel(tObs);
    x = zeros(N, 1);
    y = zeros(N, 1);

    x(1) = u(1);
    y(1) = u(2);

    % March through observation times
    for n = 2:N

        dt = tObs(n) - tObs(n-1);

        % Number of internal RK4 substeps
        nSteps = max(1, ceil(dt/hMax));
        h = dt/nSteps;

        for k = 1:nSteps
            k1 = f(u);
            k2 = f(u + 0.5*h*k1);
            k3 = f(u + 0.5*h*k2);
            k4 = f(u + h*k3);

            u = u + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        end

        x(n) = u(1);
        y(n) = u(2);
    end
end
