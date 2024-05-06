
% Describe IVP
f = @dgl;
xmin = 0;
xmax = 1;
stepsize = 0.001;
initial_conditions = [0,1];

% Perform approximation
[u, u_dot] = ForwEuler(f, stepsize, [xmin, xmax], initial_conditions);

% Get analytical solution
t = xmin:stepsize:xmax-stepsize;
y = sol(t);

% Calculate approximation error
e = y-u;

% Call draw function, see below
draw(t,u,y,e)

function y_dot = dgl(t,y)
    % Defines differential equation
    y_dot = -2 * t * (y.^2);
end

function y = sol(t)
    % Defines analytical solution
    y = 1 ./ (t.^2 + 1);
end

function draw(t,u,y,e)
    % Plot graphs
    plot(t,[u; y; abs(e)]);
    % Title
    title("Explicit-Euler")
    % Lables
    xlabel("t")
    % Legend
    legend("Approximation", "Analytical solution", "Absolute error")
end