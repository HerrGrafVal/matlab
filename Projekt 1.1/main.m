
global g
global l
g = 9.81;
l = 1;

% Stepsizes to use in approximation
h = [{1; "1 m"}, {10^-1; "$10^{-1}$ m"}, {10^-2; "$10^{-2}$ m"}, {10^-3; "$10^{-3}$ m"}, {10^-4; "$10^{-4}$ m"}, {10^-5; "$10^{-5}$ m"}];

% Create and label subplot for y(t)
func = subplot(2,1,1);
title("Mathematisches Pendel: Approximation", "Interpreter", "latex");
ylabel("$\varphi$/rad", "Interpreter", "latex");

% Create and label subplot for y_dot(t)
diff = subplot(2,1,2);
title("Mathematisches Pendel: Erste Zeitableitung", "Interpreter", "latex");
ylabel("$\dot{\varphi}$ $\cdot$ s/rad", "Interpreter", "latex");

% Finishing labels for both subplots
for i = [func, diff]
    axes(i);
    lgd = legend();
    hold(i, "on");
    set(lgd, "Interpreter", "latex");
    xlabel("$t/s$", "Interpreter", "latex");
end

for i = h 
    % Describe IVP
    Vektorfeld = @dgl;
    t0 = 0; 
    tf = 10;
    stepsize = i{1};
    phi0 = pi/2;
    phi_dot0 = 0;
    
    % Prepare parameters for function call
    tspan = t0:stepsize:tf;
    x0 = [phi0; phi_dot0];
        
    % Perform approximation
    [t, z] = SimulationEuler(Vektorfeld, tspan, x0);

    % Plot approximation
    axes(func)
    plot(t,z(1,:),"DisplayName","$\varphi$(t) Schrittweite " + i{2})

    % Plot approximation derivative
    axes(diff)
    plot(t,z(2,:),"DisplayName","$\dot{\varphi}$(t) Schrittweite " + i{2})
end

function x_dot = dgl(t,x)
    % Defines differential equation
    global g
    global l
    x_dot = [x(2); -g/l * sin(x(1))];
end
