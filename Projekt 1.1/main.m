
global g
global l
g = 9.81;
l = 1;

% Stepsizes to use in approximation
h = [{10^-2; "10^{-2}m"}, {10^-3; "10^{-3}m"}, {10^-4; "10^{-4}m"}, {10^-5; "10^{-5}m"}, {10^-6; "10^{-6}m"}];

hold on

for i = h 
    % Describe IVP
    Vektorfeld = @dgl;
    t0 = 0;
    tf = 10;
    stepsize = i{1};
    phi0 = pi/2;
    phi_dot0 = 0;
    
    tspan = t0:stepsize:tf;
    x0 = [phi0; phi_dot0];
        
    % Perform approximation
    [t, z] = SimulationEuler(Vektorfeld, tspan, x0);
    u = z(1,:);

    % Plot approximation
    plot(t,u,"DisplayName","Schrittweite " + i{2})
end

% Label plot
legend()
xlabel("t in s")
ylabel("x in m")
title("Mathematisches Pendel")

hold off

function x_dot = dgl(t,x)
    % Defines differential equation
    global g
    global l
    x_dot = [x(2); -g/l * sin(x(1))];
end
