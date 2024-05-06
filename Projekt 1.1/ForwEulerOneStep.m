function [u1, u_dot] = ForwEulerOneStep(f, h, ic)
% Returns y(t0 + h), y_dot(t0 + h)
%
% Parameters
% f = matlab function, dgl(t,y)
% h = stepsize (smallest time interval)
% ic = [t0, y0] = initial conditions
    t0 = ic(1);
    u0 = ic(:,2);
    u1 = u0 + h*f(t0, u0);
    u_dot = (u1 - u0)/h;
end