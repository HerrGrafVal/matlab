function [y, y_dot] = ForwEuler(f, h, I, ic)
% Returns y(t,y) and y_dot where y solves f
%
% Parameters
% f = matlab function, dgl(t,y)
% h = stepsize (smallest time interval)
% I = [tmin, tmax] = timeintervall
% ic = [t0, y0] = initial conditions
    tmin = I(1);
    tmax = I(2);
    t0 = ic(1);
    y0 = ic(2);

    count = (tmax - tmin) / h;
    y = zeros(1, count);
    y_dot = zeros(1, count);
    y(1) = y0;
    

    for i = 2:count
        [y(i), y_dot(i)] = ForwEulerOneStep(f, h, ic);
        ic(1) = ic(1) + h;
        ic(2) = y(i);
    end
    
end