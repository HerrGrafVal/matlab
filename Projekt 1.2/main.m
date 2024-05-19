
global g
global l
g = 9.81;
l = 1;

% Stepsizes to use in approximation
h = [{1; "1 m"}, {10^-1; "$10^{-1}$ m"}, {10^-2; "$10^{-2}$ m"}, {10^-3; "$10^{-3}$ m"}, {10^-4; "$10^{-4}$ m"}, {10^-5; "$10^{-5}$ m"}];

% Describe IVP
Vektorfeld = @dgl;
t0 = 0; 
tf = 1;
phi0 = pi/2;
phi_dot0 = 0;
global x0
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
set(e_max, "XDir", "reverse")
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
% the same calculation (as smaller stepsizes by factor 10 always include bigger ones)
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
    [tspan, z] = SimulationEuler(Vektorfeld, tspan, x0);

    % Calculate global error
    if index == 1
        % PendulumTrueSolution() doesn't work for t = [0, 1]
        % Instead only the first and last value of PendulumTrueSolution
        % with tspan = 0:0.1:1 are being used
        x_with_sufficient_stepsize = TrueSolution(0:0.1:1, x0, l, g);
        phi_with_current_stepsize = [x_with_sufficient_stepsize(1,1), x_with_sufficient_stepsize(1,end)];
        e = abs(phi_with_current_stepsize - z(1,:));
    else 
        x_with_current_stepsize = TrueSolution(tspan, x0, l, g);
        phi_with_current_stepsize = x_with_current_stepsize(1,:);
        e = abs(phi_with_current_stepsize - z(1,:));
    end    

    % Calculate local error
    tau = LocalError(tspan,stepsize,z,TrueSolution);

    % Plot global error over time
    % Since stepsize 1 corresponds to global error [0, e(1)] and
    % log10(0) = -inf this stepsize is rendered differently
    axes(e_time);
    if index == 1
        plot(tspan, log10(e), "x", "DisplayName", "e(t) Schrittweite " + i{2})
    else
        plot(tspan, log10(e), "DisplayName", "e(t) Schrittweite " + i{2})
    end

    % Note: By changing line 149 to: 
    % tau = [0];
    % And changing line 152 to:
    % for i = 1:length(t)
    % You can remove the following if/else statement and instead call:
    % plot(tspan, tau, "DisplayName", "$\tau$(t) Schrittweite " + i{2})
    % In all iterations, but since tau(0,h) = 0 for all h
    % We prefer to render local errors as follows

    % Plot local error over time
    axes(e_loc);
    if index == 1
        plot(tspan(2:end), tau, "x", "DisplayName", "$\tau$(t) Schrittweite " + i{2})
    else
        plot(tspan(2:end), tau, "DisplayName", "$\tau$(t) Schrittweite " + i{2})
    end

    % Save maximum global error
    h{3,index} = max(e);
    index = index + 1;

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

function tau = LocalError(t, h, u, sol)
% Returns tao(t_i, h) local error of approximation u for
% PendulumTrueSolution
%
% Parameters
% t = time values to evaluate
% h = stepsize
% u = approximation results
% sol = function handle to PendulumTrueSolution()

    global g
    global l
    global x0

    tau = [];
    y = sol(t, x0, l, g);

    for i = 1:length(t)-1
        % Where index 0 = t_i
        % Where index 1 = t_i + h
        t0 = t(i);
        u0 = u(1,i);

        y0 = y(:,i);
        y1 = y(1,i+1);

        f = dgl(t0, y0);
        f0 = f(2);

        tau = [tau, (y1 - u0)/h - f0];
    end


end
