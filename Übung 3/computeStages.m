function [k_mat] = computeStages(f, t_i, u_i, h, butcher)

m = length(butcher.gamma);
n_y = size(u_i, 1);

k_mat = NaN*ones(n_y, m);

for l = 1:m
    u_i_l = u_i;
    for n = 1:l-1
        u_i_l = u_i_l + h * butcher.beta(l,n) .* k_mat(:, n);
    end
    k_l = f(t_i + butcher.alpha * h, u_i_l);
    k_mat(:, l) = k_l;
end

end