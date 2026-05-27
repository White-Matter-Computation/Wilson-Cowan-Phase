% stefan ruschel, 19/05/2026
%% Extra condition delayed Wilson-Cowan restricted to sync manifold
function [res,J] = cond_dWC_sync(p, ind)

    % Hard-coded for:
    %   tau_0 = internal delay
    %   tau_1,...,tau_10 = connection delays
    %   T = period parameter

    N = 11;

    % Conditions:
    %   tau_k = min(k,N-k)*tau_1, k=2,...,10  -> 9 conditions
    %   p.period = p.parameter(ind.T)         -> 1 condition
    res = zeros(N-1,1);
    J = repmat(p_axpy(0,p,[]), N-1, 1);

    tau1 = ind.tau_1;

    for k = 2:N-1
        row = k-1;

        tauk = ind.(sprintf('tau_%d', k));
        dist_k = min(k, N-k);

        res(row) = p.parameter(tauk) - dist_k*p.parameter(tau1);

        J(row).parameter(tauk) = 1;
        J(row).parameter(tau1) = -dist_k;
    end

    % period = T
    row = N-1;

    res(row) = p.period - p.parameter(ind.T);

    J(row).period = 1;
    J(row).parameter(ind.T) = -1;
end

