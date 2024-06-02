
global g
global l
g = 9.81;
l = 1;

% Stepsizes to use in approximation
h_vec = exp(linspace(log(10^-4), log(10^1), 20));
h = cell(5, length(h_vec));
for i = 1:length(h_vec)
    h{1, i} = h_vec(i);
    h{2, i} = num2str(h_vec(i));
end

% Option to keep figure with simulation results
PLOT_APPROXIMATION = false;
STEPSIZES_TO_PLOT = [13,10,8];

% Describe IVP
Vektorfeld = @dgl;
t0 = 0; 
tf = 10;
phi0 = pi/2;
phi_dot0 = 0;
x0 = [phi0; phi_dot0];

% --------------------------------------------------------------------------

t = t0:h{1,1}:tf;
sol = PendulumTrueSolution(t, x0, l, g);

% --------------------------------------------------------------------------

close all;
% Create figure for simulation results
sim = figure();

% Create and title subplot for SimulationHeun results
heun = subplot(3,1,1);
title('Mathematisches Pendel: Simulation mit Verfahren von Heun', "Interpreter", "latex");

% Create and title subplot for SimulationCollatz results
collatz = subplot(3,1,2);
title('Mathematisches Pendel: Simulation mit Verfahren von Collatz', "Interpreter", "latex");

% Create and title subplot for SimulationRK4 results
rk4 = subplot(3,1,3);
title('Mathematisches Pendel: Simulation mit RK4 Verfahren', "Interpreter", "latex");

% Finish labeling
for ax = [heun, collatz, rk4]
    axes(ax);
    xlabel("t in s", "Interpreter", "latex");
    ylabel("x in rad", "Interpreter", "latex");
    lgd = legend();
    set(lgd, "Interpreter", "latex");
    hold(ax, "on");
    plot(t, sol(1,:), "DisplayName", 'Analytische L\"osung')
end

% --------------------------------------------------------------------------

% Create figure for method errors
err = figure();

% Create and title subplot for SimulationHeun error
err_heun = subplot(3,1,1);
title('Maximaler Fehler pro Schrittweite f\"ur Verfahren von Heun', "Interpreter", "latex");
ylabel("$max_{t}|e_{Heun}(h)|$", "Interpreter", "latex");

% Create and title subplot for SimulationCollatz error
err_collatz = subplot(3,1,2);
title('Maximaler Fehler pro Schrittweite f\"ur Verfahren von Collatz', "Interpreter", "latex");
ylabel("$max_{t}|e_{Collatz}(h)|$", "Interpreter", "latex");

% Create and title subplot for SimulationHeun error
err_rk4 = subplot(3,1,3);
title('Maximaler Fehler pro Schrittweite f\"ur RK4 Verfahren', "Interpreter", "latex");
ylabel("$max_{t}|e_{RK4}(h)|$", "Interpreter", "latex");

% Finish labeling
for ax = [err_heun, err_collatz, err_rk4]
    axes(ax);
    xlabel("$log_{10}(h)$", "Interpreter", "latex");
    lgd = legend();
    set(lgd, "Interpreter", "latex");
    hold(ax, "on");
end

% --------------------------------------------------------------------------

methods = [{@SimulationHeun; heun; err_heun}, {@SimulationCollatz; collatz; err_collatz}, {@SimulationRK4; rk4; err_rk4}];

% --------------------------------------------------------------------------

% Repeat for all stepsizes
figure(sim);
index = 1;
for i = h 
    stepsize = i{1};

    % Prepare parameters for function call
    tspan = t0:stepsize:tf;

    % Calculate true solution with current stepsize
    x_with_current_stepsize = PendulumTrueSolution(tspan, x0, l, g);
    phi_with_current_stepsize = x_with_current_stepsize(1,:);
    
    subindex = 3;
    for method = methods
        % Perform approximation
        [tspan, z] = method{1}(Vektorfeld, tspan, x0);

        % Plot approximation result for desired stepsizes
        if ismember(index, STEPSIZES_TO_PLOT)
            axes(method{2});
            plot(tspan, z(1,:), "DisplayName", "Schrittweite " + i{2});
        end

        % Calculate and save error maximum
        e = abs(phi_with_current_stepsize - z(1,:));
        h{subindex, index} = max(e);
        subindex = subindex + 1;
    end

    index = index + 1;
end

% --------------------------------------------------------------------------

% Repeat for all methods
figure(err);
stepsizes = cell2mat(h(1, :));
log10_h = log10(stepsizes);
index = 3;
for method = methods
    % Plot maximum error over stepsize
    e = cell2mat(h(index, :));
    log10_e = log10(e);
    axes(method{3});
    plot(log10_h, log10_e, "x", "DisplayName", "Maximaler Fehler")
    index = index + 1;

    % Interpolate linear function for log10_e over log10_h
    e_of_h_coefficients = polyfit(log10_h, log10_e, 1);
    h_ = [-4, 1];
    e_ = polyval(e_of_h_coefficients, h_);
    plot(h_, e_, "DisplayName", "Lineare Regressionsgerade")

    % Annotate linear function
    text = "$log_{10}(e)$ = " + num2str(e_of_h_coefficients(1)) + " $\cdot$ $log_{10}(h)$ + " + num2str(e_of_h_coefficients(2) + " ");
    arrow_x = [-1.75, -1.5];
    arrow_y = [0, polyval(e_of_h_coefficients, -1.5)]
    anno = annotation("textarrow", "String", text, "HeadStyle", "none", "Interpreter", "latex");
    anno.Parent = method{3};
    anno.X = arrow_x;
    anno.Y = arrow_y;
end

% --------------------------------------------------------------------------

% Close figure containing simulation results if PLOT_APPROXIMATION is set to false
if PLOT_APPROXIMATION == false
    close(sim)
end

% --------------------------------------------------------------------------

function x_dot = dgl(t,x)
    % Defines differential equation
    global g
    global l
    x_dot = [x(2); -g/l * sin(x(1))];
end
