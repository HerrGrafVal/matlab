
% Stepsizes to use in approximation
h_vec = exp(linspace(log(10^-4), log(10^0), 10));
l = length(h_vec);
h = cell(2,l);
for i = 1:l
    h(1,i) = {h_vec(i)};
    h(2,i) = {num2str(h_vec(i))};
end

% Select which stepsizes to plot
start = 7;
finish = 10;

% Label plot for y(t), u(t)
name = ["Verfahren von Heun"; "$\dot{y}(t) = y(t) + t$"];
title(name, "Interpreter", "latex");
ylabel("y(t)", "Interpreter", "latex");
hold on;
lgd = legend();
set(lgd, "Interpreter", "latex");
xlabel("t", "Interpreter", "latex");

for i = h(:,start:finish)
    % Describe IVP
    Vektorfeld = @dgl;
    ts = 0; 
    tf = 8;
    stepsize = i{1};
    
    % Prepare parameters for function call
    t0 = 1;
    y0 = 2;
    tspanb = ts:stepsize:t0;
    tspanf = t0:stepsize:tf;
        
    % Perform approximation
    [t_b, u_b] = SimulationHeunBackward(Vektorfeld, tspanb, y0);
    [t_f, u_f] = SimulationHeunForward(Vektorfeld, tspanf, y0);

    % Plot approximation
    plot([t_b, t_f], [u_b, u_f],"DisplayName","u(t) Schrittweite " + i{2})
end

% Plot analytical solution
t = ts:h_vec(1):tf;
plot(t, sol(t), "DisplayName", 'y(t) Analytische L\"osung')

function y_dot = dgl(t,y)
    % Defines differential equation
    y_dot = y + t;
end

function y = sol(t)
    % Defines analytical solution
    y = 4*exp(t-1) -t -1;
end
