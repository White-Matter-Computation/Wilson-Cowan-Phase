clc; clear all; close all;
addpath('../DDE-Biftool/ddebiftool',...
    '../DDE-Biftool/ddebiftool_extra_psol',...
    '../DDE-Biftool/ddebiftool_extra_nmfm',...
    '../DDE-Biftool/ddebiftool_extra_symbolic',...
    '../DDE-Biftool/ddebiftool_utilities');

data   = load('data/dde_batch_0.2_0.0_eps_0.1_20.mat');
omegas = data.omegas;
tau1s  = data.tau1s;
n_runs = numel(omegas);

results = cell(n_runs, 1);
tol_zero = 1e-6;   % |Re(λ)| < tol → on boundary
status   = zeros(n_runs, 1);
num_uns = zeros(n_runs,1);
for i = 1:n_runs
    A   = data.A_coeff{i};   % N x N x n_delays_i
    Tau = data.Tau{i};        % [0, tau_0, theta_1, ...]

    lam = dde_stst_eig_cheb(A, Tau, ...
        'minimal_real_part',         -2,   ...
        'root_accuracy',           1e-8,   ...
        'max_number_of_eigenvalues',  30,  ...
        'min_number_of_eigenvalues',  20);

    results{i} = lam.l0;
    re = real(lam.l0);

    
    %can make 0 to see all have the trivial 0 eigenvalue
    if sum(abs(re) < tol_zero) > 1
        status(i) = 2;          % bif
    elseif any(re > tol_zero)
        status(i) = 1;          % unstable
        num_uns(i) = sum(re>tol_zero);
    else
        status(i) = 0;          % stable
    end
end

%% Plot all spectra coloured by omega
figure(1); clf;

cmap = jet(256);

% Precompute ranges (avoids recomputing + avoids divide-by-zero)
omega_min = min(omegas);
omega_range = max(omegas) - omega_min + eps;

tau_min = min(tau1s);
tau_range = max(tau1s) - tau_min + eps;

% Plot 1: coloured by Omega
subplot(1,2,1); hold on;
for i =1:n_runs
    c = (omegas(i) - omega_min) / omega_range;
    idx = 1 + round(c * 255);
    plot(real(results{i}), imag(results{i}), '.', 'Color', cmap(idx, :));
end
colormap(cmap);
clim([omega_min omega_min + omega_range]);
cb = colorbar;
ylabel(cb, '\Omega');
xline(0,'--r'); xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
title('Coloured by \Omega');

% Plot 2: coloured by tau_1
subplot(1,2,2); hold on;
for i = 1:n_runs
    c = (tau1s(i) - tau_min) / tau_range;
    idx = 1 + round(c * 255);
    plot(real(results{i}), imag(results{i}), '.', 'Color', cmap(idx, :));
end
colormap(cmap);
clim([tau_min tau_min + tau_range]);
cb = colorbar;
ylabel(cb, '\tau_1');
xline(0,'--r'); xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
title('Coloured by \tau_1');

% 
% %% Split by status
idx_stable   = status == 0;
idx_unstable = status == 1;
idx_boundary = status == 2;

%% Plot Omega vs tau_1
figure(2); clf; hold on;

scatter(tau1s(idx_stable),omegas(idx_stable),     40, ...
    'o', ...
    'MarkerFaceColor', [0.2 0.6 1.0], ...
    'MarkerEdgeColor', 'none', ...
    'DisplayName', 'Stable');

scatter( tau1s(idx_unstable),omegas(idx_unstable), 40, ...
    's', ...
    'MarkerFaceColor', [1.0 0.3 0.3], ...
    'MarkerEdgeColor', 'none', ...
    'DisplayName', 'Unstable');

scatter( tau1s(idx_boundary),omegas(idx_boundary), 80, ...
    'd', ...
    'MarkerFaceColor', [0.6 0.0 0.8], ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Boundary (Re(\lambda)\approx 0)');

xlabel('\tau_1');
ylabel('\Omega');
title('Stability map in (\Omega, \tau_1) plane');
legend('Location', 'best');
grid on; box on;
hold off;

%% Optional: print summary
fprintf('Stable:   %d / %d\n', sum(idx_stable),   n_runs);
fprintf('Unstable: %d / %d\n', sum(idx_unstable), n_runs);
fprintf('Boundary: %d / %d\n', sum(idx_boundary), n_runs);

figure(3); clf;
plot(real(results{n_runs}), imag(results{n_runs}), 'ob', 'MarkerFaceColor', 'b');
xline(0, '--r'); yline(0, '--k');
xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
title(sprintf('Spectrum (N=%d, %d delays)', size(A,1), numel(Tau)-1));
grid on;



%%
num_map = unique(num_uns);
cmap = lines(length(num_map));
figure(4); clf;
hold on

for i = 1:length(num_map)
    j = num_map(i);
    idx = (j == num_uns);
    scatter(tau1s(idx), omegas(idx), 10, cmap(i,:), 'filled');
end

hold off

colormap(cmap);
cb = colorbar;
cb.Ticks = 1:length(num_map);
cb.TickLabels = num_map;   % show actual j values

xlabel('Re(\lambda)');
ylabel('Im(\lambda)');
xline(0,'--r');
title('Coloured by \tau_1');

