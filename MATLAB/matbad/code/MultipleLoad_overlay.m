clc; clear all; close all;
addpath('../DDE-Biftool/ddebiftool',...
    '../DDE-Biftool/ddebiftool_extra_psol',...
    '../DDE-Biftool/ddebiftool_extra_nmfm',...
    '../DDE-Biftool/ddebiftool_extra_symbolic',...
    '../DDE-Biftool/ddebiftool_utilities');

%% Define datasets 
dataset_files  = {'data/dde_batch_0.2_0.0_eps_0.01.mat', ...
                  'data/dde_batch_0.2_0.57_eps_0.01.mat', ...
                  'data/dde_batch_0.2_1.14_eps_0.01.mat', ...
                  'data/dde_batch_0.2_1.71_eps_0.01.mat', ...
                  'data/dde_batch_0.2_2.28_eps_0.01.mat', ...
                  'data/dde_batch_0.2_2.86_eps_0.01.mat'};



dataset_labels  = {'\psi=0.00', '\psi=0.57', '\psi=1.14', ...
                   '\psi=1.71', '\psi=2.28', '\psi=2.86'};
dataset_markers = {'o', 's', '^', 'd', 'v', 'p'};   % 6 distinct shapes

status_colors = [0.2 0.6 1.0;    % stable   → blue
                 1.0 0.3 0.3;    % unstable → red
                 0.6 0.0 0.8];   % boundary → purple

%% Encoding — line style per status, color per dataset
dataset_colors = [0.2 0.4 0.9;    % psi=0.00  → blue
                  0.9 0.4 0.1;    % psi=0.57  → orange
                  0.1 0.7 0.3;    % psi=1.14  → green
                  0.8 0.1 0.5;    % psi=1.71  → magenta
                  0.1 0.7 0.8;    % psi=2.28  → cyan
                  0.5 0.3 0.8];   % psi=2.86  → purple

status_styles = {'-', '--', ':'};   % stable, unstable, boundary
status_widths = [2.0, 2.0, 1.5];
status_labels = {'Stable', 'Unstable', 'Boundary'};
tol_zero      = 1e-6;

%% Verify files eixst 


%% Process each dataset (unchanged)
n_datasets = numel(dataset_files);
all_status = cell(n_datasets, 1);
all_omegas = cell(n_datasets, 1);
all_tau1s  = cell(n_datasets, 1);

% for d = 1:n_datasets
%     if exist(dataset_files{d}, 'file')
%         fprintf('OK:      %s\n', dataset_files{d});
%     else
%         fprintf('MISSING: %s\n', dataset_files{d});
%     end
% end


for d = 1:n_datasets
    data   = load(dataset_files{d});
    omegas = data.omegas(:);
    tau1s  = data.tau1s(:);
    n_runs = numel(omegas);
    status = zeros(n_runs, 1);

    for i = 1:n_runs
        A   = data.A_coeff{i};
        Tau = data.Tau{i};
        lam = dde_stst_eig_cheb(A, Tau, ...
            'minimal_real_part',          -2,  ...
            'root_accuracy',            1e-8,  ...
            'max_number_of_eigenvalues',   30, ...
            'min_number_of_eigenvalues',   20);
        re = real(lam.l0);
        if sum(abs(re) < tol_zero) > 1
            status(i) = 2;
        elseif any(re > tol_zero)
            status(i) = 1;
        end
    end

    all_status{d} = status;
    all_omegas{d} = omegas;
    all_tau1s{d}  = tau1s;
    fprintf('Loaded %s: %d runs\n', dataset_files{d}, n_runs);
end

%% Plot
%% Plot
figure(1); clf;
set(gcf, 'Color', 'w', ...          % white figure background
         'Position', [100 100 1100 600]);   % larger figure
hold on;

% White axes background + thicker grid
ax = gca;
ax = gca;
ax.Color            = 'w';
ax.XColor           = 'k';   % x-axis line + tick labels
ax.YColor           = 'k';   % y-axis line + tick labels
ax.GridColor        = [0.7 0.7 0.7];
ax.GridAlpha        = 0.8;
ax.GridLineStyle    = '-';
ax.LineWidth        = 1.2;
ax.FontSize         = 13;
ax.TickDir          = 'out';
ax.Box              = 'on';

% Bolder dataset colors — higher contrast against white
dataset_colors = [0.0 0.3 0.8;    % psi=0.00  → deep blue
                  0.9 0.4 0.0;    % psi=0.57  → deep orange
                  0.1 0.6 0.1;    % psi=1.14  → deep green
                  0.8 0.0 0.4;    % psi=1.71  → deep pink
                  0.0 0.6 0.7;    % psi=2.28  → teal
                  0.5 0.1 0.7];   % psi=2.86  → deep purple

% Thicker lines + more distinct styles
status_styles = {'-', '--', ':'};
status_widths = [1.8, 1.5, 1.5];

% Seed dataset legend
for d = 1:n_datasets
    h_dataset(d) = plot(nan, nan, '-', ...
        'Color',       dataset_colors(d,:), ...
        'LineWidth',   1.8, ...
        'DisplayName', dataset_labels{d});
end

% Seed status legend
for s = 0:2
    h_status(s+1) = plot(nan, nan, status_styles{s+1}, ...
        'Color',       [0.1 0.1 0.1], ...
        'LineWidth',   status_widths(s+1), ...
        'DisplayName', status_labels{s+1});
end

% Plot segments
for d = 1:n_datasets
    omegas = all_omegas{d};
    tau1s  = all_tau1s{d};
    status = all_status{d};
    n_pts  = numel(omegas);
    col    = dataset_colors(d,:);

    i = 1;
    while i <= n_pts
        s       = status(i);
        j       = i;
        while j <= n_pts && status(j) == s
            j = j + 1;
        end
        seg_end = min(j, n_pts);
        plot(tau1s(i:seg_end), omegas(i:seg_end), ...
            status_styles{s+1}, ...
            'Color',            col, ...
            'LineWidth',        status_widths(s+1), ...
            'HandleVisibility', 'off');
        i = j;
    end
end

% Legend outside the plot so it doesn't obscure data
lg = legend([h_dataset(:); h_status(:)], ...
    [dataset_labels, status_labels], ...
    'Location', 'eastoutside', ...
    'NumColumns', 1, ...
    'FontSize', 12, ...
    'TextColor', 'k', ...
    'Color', 'w', ...
    'Box', 'on');
lg.BoxFace.ColorType = 'truecoloralpha';
lg.BoxFace.ColorData = uint8([255 255 255 230]');   % semi-transparent white

% xlabel('\tau_1',  'FontSize', 14, 'FontWeight', 'bold');
% ylabel('\Omega',  'FontSize', 14, 'FontWeight', 'bold');
% title('Stability map in (\Omega, \tau_1) plane', 'FontSize', 14);
xlabel('\tau_1', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
ylabel('\Omega', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
title('Stability map in $(\Omega, \tau_1)$ plane for $\tau_0=0.2$ and $\varepsilon=0.001$', ...
      'FontSize', 14, ...
      'Color', 'k', ...
      'Interpreter', 'latex');
ax.XTick = 0:1:10;
grid on; box on; hold off;

% Export as high-res PNG
exportgraphics(gcf, 'stability_map_0_01.png', 'Resolution', 300, 'BackgroundColor', 'white');


%%
%%Get the start and end points of that 2nd stable branch 
%% Get transition points for dataset 1
status_1 = all_status{1};
omegas_1 = all_omegas{1};
tau1s_1  = all_tau1s{1};

for i = 1:length(status_1)-1
    if status_1(i) == 0 && status_1(i+1) == 1
        fprintf('Stable → Unstable at index %d:  tau1=%.4f  Omega=%.4f\n', ...
            i, tau1s_1(i), omegas_1(i));
    elseif status_1(i) == 1 && status_1(i+1) == 0
        fprintf('Unstable → Stable at index %d:  tau1=%.4f  Omega=%.4f\n', ...
            i+1, tau1s_1(i+1), omegas_1(i+1));
    end
end