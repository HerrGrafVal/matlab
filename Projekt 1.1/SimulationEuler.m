function [t, z] = SimulationEuler(Vektorfeld, tspan, x0)
% Returns y(t,y) and y_dot where y solves f
%
% Parameters
% Vektorfeld = matlab function, dgl(t,y)
% tspan = t0:stepsize:tf = timeintervall
% x0 = initial condition at t0

    h = tspan(2)-tspan(1);
    t0 = tspan(1);

    u = zeros(size(x0, 1), length(tspan));
    u_dot = u;
    u(:,1) = x0;
    ic = [[t0; t0], x0];

    for i = 2:length(tspan)
        [u(:,i), u_dot(:,i)] = ForwEulerOneStep(Vektorfeld, h, ic);
        ic(:,1) = ic(:,1) + h;
        ic(:,2) = u(:,i);
    end

    t = tspan;
    z = u;
    
end