function [t, u] = SimulationHeunBackward(Vektorfeld, tspan, y0)
% Returns [t, z] where z(t) solves Vektorfeld
%
% Parameters
% Vektorfeld = matlab function, dgl(t,y)
% tspan = tstart:stepsize:t0 = timeintervall
% y0 = initial condition at t0

    % Flip time vector for backwards iteration
    t = tspan;
    tspan = flip(tspan);

    % Setup
    h = tspan(2)-tspan(1);
    t0 = tspan(1);
    u0 = y0;
    f = Vektorfeld;

    u = zeros(size(y0, 1), length(tspan));
    u(:,1) = u0;

    for i = 2:length(tspan)
        % Approximate
        t1 = t0 + h;
        l1 = u0 + h*f(t0, u0);
        u(:,i) = u0 + 0.5 * h * (f(t0, u0) + f(t1, l1));

        % Setup for next iteration
        t0 = t0 + h;
        u0 = u(:,i);
    end

    % Return to chronological order
    u = flip(u);
    
end