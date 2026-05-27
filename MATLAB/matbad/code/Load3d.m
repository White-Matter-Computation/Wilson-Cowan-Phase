% Spectrum of linearised problem - multi tau_0 version

clc; clear all; close all; 
% Loading the path where the DDE-Biftool scripts are stored
addpath('../DDE-Biftool/ddebiftool',...
    '../DDE-Biftool/ddebiftool_extra_psol',...
    '../DDE-Biftool/ddebiftool_extra_nmfm',...
    '../DDE-Biftool/ddebiftool_extra_symbolic',...
    '../DDE-Biftool/ddebiftool_utilities');

data   = load('dde_batch_3d.mat');
tau0s  = data.tau0s(:);
n_tau0 = numel(tau0s);
tol_zero = 1e-6;

%% Compute stability for all (tau_0, run) pairs
% Store as cells since each branch has different length
status_all = cell(n_tau0, 1);
omegas_all = cell(n_tau0, 1);
tau1s_all  = cell(n_tau0, 1);

for k = 1:n_tau0
    omegas_k = data.omegas{k}(:);
    tau1s_k  = data.tau1s{k}(:);
    n_runs   = numel(omegas_k);
    status_k = zeros(n_runs, 1);

    for i = 1:n_runs
        A   = data.A_coeff{k}{i};
        Tau = data.Tau{k}{i};
        lam = dde_stst_eig_cheb(A, Tau, ...
            'minimal_real_part',          -2,  ...
            'root_accuracy',            1e-8,  ...
            'max_number_of_eigenvalues',   30, ...
            'min_number_of_eigenvalues',   20);
        re = real(lam.l0);
        if sum(abs(re) < tol_zero) > 1
            status_k(i) = 2;       % bifurcation boundary
        elseif any(re > tol_zero)
            status_k(i) = 1;       % unstable
        else
            status_k(i) = 0;       % stable
        end
    end

    status_all{k} = status_k;
    omegas_all{k} = omegas_k;
    tau1s_all{k}  = tau1s_k;
    fprintf('tau_0=%.4f done (%d runs)\n', tau0s(k), n_runs);
end

%% Colours and markers
colors = [0.2 0.6 1.0;   % stable   → blue
          1.0 0.3 0.3;   % unstable → red
          0.6 0.0 0.8];  % boundary → purple
labels  = {'Stable', 'Unstable', 'Boundary'};
markers = {'o', 's', 'd'};
sizes   = [30, 30, 60];

%% Figure 1: 3D scatter (Omega, tau_1, tau_0)
%% Figure 1: 3D scatter — Omega on Z axis
figure(1); clf; hold on;
legend_added = false(3,1);
% Check data is sensible
fprintf('n_tau0 = %d\n', n_tau0);
for k = 1:n_tau0
    fprintf('k=%d: %d points, status counts: stable=%d unstable=%d boundary=%d\n', ...
        k, numel(omegas_all{k}), ...
        sum(status_all{k}==0), ...
        sum(status_all{k}==1), ...
        sum(status_all{k}==2));
end

for k = 1:n_tau0
    omegas_k = omegas_all{k};
    tau1s_k  = tau1s_all{k};
    status_k = status_all{k};
    tau0_k   = tau0s(k) * ones(size(omegas_k));

    for s = 0:2
        idx = status_k == s;
        if ~any(idx); continue; end
        sc = scatter3(tau1s_k(idx), tau0_k(idx), omegas_k(idx), ...  % ← reordered
            sizes(s+1), ...
            'Marker',          markers{s+1}, ...
            'MarkerFaceColor', colors(s+1,:), ...
            'MarkerEdgeColor', 'none', ...
            'DisplayName',     labels{s+1});
        if legend_added(s+1)
            sc.HandleVisibility = 'off';
        else
            legend_added(s+1) = true;
        end
    end
end

xlabel('\tau_1'); ylabel('\tau_0'); zlabel('\Omega');  % ← updated
title('Stability in (\tau_1, \tau_0, \Omega) space');
legend('Location', 'best');
grid on; box on;
view(45, 25);
hold off;
%% Figure 2: 2D top-down view — tau_1 vs tau_0 only
% For each (tau_0, tau_1) point take the worst status across all Omega values
figure(2); clf; hold on;
legend_added2 = false(3,1);

for k = 1:n_tau0
    tau1s_k  = tau1s_all{k};
    status_k = status_all{k};

    % Unique tau_1 values on this branch, take max status (worst case)
    tau1_unique = unique(tau1s_k);
    for m = 1:numel(tau1_unique)
        mask        = tau1s_k == tau1_unique(m);
        worst_status = max(status_k(mask));   % 0<1<2, boundary dominates
        sc = scatter(tau1_unique(m), tau0s(k), ...
            sizes(worst_status+1), ...
            'Marker',          markers{worst_status+1}, ...
            'MarkerFaceColor', colors(worst_status+1,:), ...
            'MarkerEdgeColor', 'none', ...
            'DisplayName',     labels{worst_status+1});
        if legend_added2(worst_status+1)
            sc.HandleVisibility = 'off';
        else
            legend_added2(worst_status+1) = true;
        end
    end
end

xlabel('\tau_1'); ylabel('\tau_0');
title('Stability map in (\tau_1, \tau_0) plane');
legend('Location', 'best');
grid on; box on;
hold off;