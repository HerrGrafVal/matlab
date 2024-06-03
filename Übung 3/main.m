butcher.alpha = [0; 1/2; 1];
butcher.beta = [0 0 0; 1/2 0 0; -1 2 0];
butcher.gamma = [1/6 4/6 1/6];

function y_dot = dgl(t, y)
    y_dot = -y;
end