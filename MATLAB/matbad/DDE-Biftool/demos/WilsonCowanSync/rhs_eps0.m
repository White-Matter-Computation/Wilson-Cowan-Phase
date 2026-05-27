function dy = rhs_eps0(~, y, Z, p, ind)
    % Current state
    u = y(1);
    v = y(2);

    % Delayed state at t - tau_0
    u_tau0 = Z(1,1);
    v_tau0 = Z(2,1);

    z_u = p(ind.I_u) ...
        + p(ind.w_ee)*u_tau0 ...
        - p(ind.w_ei)*v_tau0;

    z_v = p(ind.I_v) ...
        + p(ind.w_ie)*u_tau0 ...
        - p(ind.w_ii)*v_tau0;

    dy = zeros(2,1);
    beta = 20;
    S = @(z,beta) 1 ./ (1 + exp(-beta*z));
    dy(1) = (-u + S(z_u, beta)) / p(ind.tau_e);
    dy(2) = (-v + S(z_v, beta)) / p(ind.tau_i);
end
