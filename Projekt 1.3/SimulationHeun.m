function [t, z] = SimulationHeun(Vektorfeld, tspan, x0)
% Returns [t, z] where z(t) solves Vektorfeld
%
% Parameters
% Vektorfeld = matlab function handle, dgl(t,y)
% tspan = t0:stepsize:tend = timeintervall
% y0 = initial condition at tspan(1)
   
    % Setup
    t = tspan;

    h = tspan(2)-tspan(1);
    t0 = tspan(1);
    z0 = x0;
    f = Vektorfeld;

    z = zeros(size(x0,1), length(tspan));
    z(:,1) = z0;

    for i = 2:length(t)
        % Approximate according to Heun's method
        k1 = f(t0, z0);
        k2 = f(t0 + h, z0 + h * k1);
        z(:,i) = z0 + h * (0.5 * k1 + 0.5 * k2);

        % Setup for next iteration 
        t0 = t0 + h;
        z0 = z(:,i);
    end
    
end