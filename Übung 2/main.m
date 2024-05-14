
% Stepsizes to use in approximation
h_vec = exp(linspace(log(10^-4), log(10^0), 10));
h = [h_vec; zeros(size(h_vec))];
index = 1;

% Label plot for e(h)
name = ["Verfahren von Heun"; "$\dot{y}(t) = y(t) + t$"];
title(name, "Interpreter", "latex");
ylabel("e(h)", "Interpreter", "latex");
hold on;
lgd = legend();
set(lgd, "Interpreter", "latex");
xlabel("$ln(h)$", "Interpreter", "latex");

for i = h(1,:)
    % Describe IVP
    Vektorfeld = @dgl;
    ts = 0; 
    tf = 8;
    stepsize = i;
    
    % Prepare parameters for function call
    t0 = 1;
    y0 = 2;
    tspanb = ts:stepsize:t0;
    tspanf = t0:stepsize:tf;
        
    % Perform approximation
    [t_b, u_b] = SimulationHeunBackward(Vektorfeld, tspanb, y0);
    [t_f, u_f] = SimulationHeunForward(Vektorfeld, tspanf, y0);
    t = [t_b, t_f];
    u = [u_b, u_f];

    % Calculate error
    y = sol(t);
    e = abs(u - y);

    % Save maximum error per stepwidth
    h(2, index) = max(e);
    index = index + 1;
end

% Plot e(h)
log_h = log(h(1,:));
plot(log_h, h(2,:), "xr", "DisplayName", 'e(h) Maximum error per stepwidth')

function y_dot = dgl(t,y)
    % Defines differential equation
    y_dot = y + t;
end

function y = sol(t)
    % Defines analytical solution
    y = 4*exp(t-1) -t -1;
end
