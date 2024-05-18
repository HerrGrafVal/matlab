function u1 = ForwEulerOneStep(f, h, t0, u0)
% Returns u1
%
% Parameters
% f = matlab function, dgl(t,y)
% h = stepsize (smallest time interval)
% t0, u0 = initial conditions
    u1 = u0 + h*f(t0, u0);
    %u_dot = (u1 - u0)/h;
end