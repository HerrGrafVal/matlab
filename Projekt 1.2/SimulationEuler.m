function [t, z] = SimulationEuler(Vektorfeld, tspan, x0)
% Returns [t, z] where z(t) solves Vektorfeld
%
% Parameters
% Vektorfeld = matlab function handle, dgl(t,y)
% tspan = t0:stepsize:tf = timeintervall
% x0 = initial condition at t0

    h = tspan(2)-tspan(1);
    t0 = tspan(1);

    z = zeros(size(x0, 1), length(tspan));
    z(:,1) = x0;

    for i = 2:length(tspan)

        % Display progress in %
        % if mod(i,100000) == 0
        %     disp(i/length(tspan) * 100);
        % end

        z(:,i) = ForwEulerOneStep(Vektorfeld, h, t0, x0);
        t0 = t0 + h;
        x0 = z(:,i);
    end

    t = tspan;
    
end