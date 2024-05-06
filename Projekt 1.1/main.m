
%h = [10^-2, 10^-4, 10^-8, 10^-12, 10^-16, 10^-18];
h = [10^-2, 10^-4];

global g
g = 9.81;
global l
l = 1;
hold on

for i = h  
    % Describe IVP
    Vektorfeld = @dgl;
    t0 = 0;
    tf = 10;
    stepsize = i;
    phi0 = pi/2;
    phi_dot0 = 0;
    
    tspan = t0:stepsize:tf;
    x0 = [phi0; phi_dot0];
        
    % Perform approximation
    [t, z] = SimulationEuler(Vektorfeld, tspan, x0);
    u = z(1,:);
    draw(t,u)
end


% Call draw function, see below
draw(tspan,u)

function x_dot = dgl(t,x)
    % Defines differential equation
    global g
    global l
    x_dot = [0; 0];
    x_dot(1) = x(2);
    x_dot(2) = -g/l * sin(x(1));
end

function draw(t,u)
    % Plot graphs
    plot(t, u);
    % Title
    title("Mathematisches Pendel")
    % Lables
    xlabel("t in s")
    ylabel("x in m")
    % Legend
    %legend(int2str(res))
end