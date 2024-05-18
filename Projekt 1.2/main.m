
global g
global l
g = 9.81;
l = 1;

% Stepsizes to use in approximation
h = [{1; "1 m"}, {10^-1; "$10^{-1}$ m"}, {10^-2; "$10^{-2}$ m"}, {10^-3; "$10^{-3}$ m"}];

% Describe IVP
Vektorfeld = @dgl;
t0 = 0; 
tf = 1;
phi0 = pi/2;
phi_dot0 = 0;
x0 = [phi0; phi_dot0];

% --------------------------------------------------------------------------

figure()

% Create and label subplot for global error e(t)
e_time = subplot(3,1,1);
title('Mathematisches Pendel: Globaler Fehler \"uber der Zeit', "Interpreter", "latex");
xlabel("$t/s$", "Interpreter", "latex");
ylabel("$log_{10}(e(t))$", "Interpreter", "latex");
lgd = legend();
set(lgd, "Interpreter", "latex");
hold(e_time, "on");

% Create and label subplot for global error maximum
e_max = subplot(3,1,2);
title("Maxima des globalen Fehlers pro Schrittweite", "Interpreter", "latex");
xlabel("$log_{10}(h)$", "Interpreter", "latex");
ylabel("$max_t|e(t)|$", "Interpreter", "latex");
hold(e_max, "on");

% Create and label subplot for local error 
e_loc = subplot(3,1,3);
title('Lokaler Fehler \"uber der Zeit', "Interpreter", "latex");
xlabel("$t/s$", "Interpreter", "latex");
ylabel("$\tau$(t)", "Interpreter", "latex");
lgd = legend();
set(lgd, "Interpreter", "latex");
hold(e_loc, "on");


% --------------------------------------------------------------------------

% Using cache for PendulumTrueSolution to prevent multiple executions of
% the same calculation
global TrueSolution
TrueSolution = memoize(@PendulumTrueSolution);

% --------------------------------------------------------------------------

% Repeat for all stepsizes
index = 1;
for i = h 
    stepsize = i{1};

    % Prepare parameters for function call
    tspan = t0:stepsize:tf;
    
    % Perform approximation
    % Write output to [, z] since tspan is known
    [tspan, z] = SimulationEuler(Vektorfeld, tspan, x0);

    % Calculate global error 
    x_with_current_stepsize = TrueSolution(tspan, x0, l, g);
    phi_with_current_stepsize = x_with_current_stepsize(1,:);
    e = abs(phi_with_current_stepsize - z(1,:));

    % Calculate local error
    tau = LocalError(tspan,stepsize,z(1,:));

    % Plot global error over time
    axes(e_time);
    plot(tspan, log10(e), "DisplayName", "e(t) Schrittweite " + i{2})

    % Save maximum global error
    h{3,index} = max(e);
    index = index + 1;

    % Plot local error over time
    axes(e_loc);
    plot(tspan, tau, "DisplayName", "$\tau$(t) Schrittweite " + i{2})
end

% --------------------------------------------------------------------------

axes(e_max);
plot(log10([h{1,:}]), [h{3,:}], "x");
legend("Maximaler Fehler pro Schrittweite h", "Interpreter", "latex");

% --------------------------------------------------------------------------

function x_dot = dgl(t,x)
    % Defines differential equation
    global g
    global l
    x_dot = [x(2); -g/l * sin(x(1))];
end

function [x] = PendulumTrueSolution(t, x0, l , g)
    % Analytical solution for dgl(t,x) above
    % Source: https://de.wikipedia.org/wiki/Mathematisches_Pendel
    x(1, :) = 2 * atan(tan(x0(1) / 2) * jacobiCN(sqrt(g/l) * t, sin(x0(1) / 2)));
    % x(2, :) = -1 * sqrt(g/l) * sin(x0(1)) * jacobiSD(sqrt(g/l) * t, sin(x0(1) / 2));
end

function tau = LocalError(t0, h, u_t0)
% Returns tao(t_i, h) local error of approximation u for
% PendulumTrueSolution
%
% Parameters
% t0 = t_i
% h = stepsize
% u_t0 = approximation value u_i
    global g
    global l
    global TrueSolution

    next = TrueSolution(t0 + h, u_t0, l, g);
    % First row of next contains y(t0 + h)
    y_t1 = next(1,:);

    increment = dgl(t0, TrueSolution(t0, u_t0, l, g));
    % Second row of increment contains f(t0,y(t0))
    f_t0 = increment(2,:);

    tau = (y_t1 - u_t0)/h - f_t0;
end
