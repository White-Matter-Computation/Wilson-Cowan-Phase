%% Delayed WC network resctricted to sync manifold
%
% tau_e du/dt = -u + f(I_u + w_ee u(t-tau_0)
%                          - w_ei v(t-tau_0)
%                          + eps_c sum_{j=1}^{N-1} w_j u(t-tau_j))
%
% tau_i dv/dt = -v + f(I_v + w_ie u_i(t-tau_0)
%                          - w_ii v_i(t-tau_0))
%
% f(x) = 1/(1+exp(-beta*x))
%
%% load DDE-Biftool into path
clear
base=[pwd(),'/../../'];
addpath([base,'ddebiftool/'],...
    [base,'ddebiftool_extra_psol/'],...
    [base,'ddebiftool_utilities/']);
format compact
%% Initial parameters and state
N=11;
parnames = { ...
    'tau_e', ...
    'tau_i', ...
    'I_u', ...
    'I_v', ...
    'w_ee', ...
    'w_ei', ...
    'w_ie', ...
    'w_ii', ...
    'eps_c', ...
    'T',...      % placeholder for period
    'tau_0'};

for k = 1:N-1
    parnames{end+1} = sprintf('tau_%d', k); 
end

cind = [parnames; num2cell(1:length(parnames))];
ind = struct(cind{:}); 

% initial parameters (interaction delays set to 0.0)
w_ee = 1.0;
w_ei = 2.0;
w_ie = 1.0;
w_ii = 0.25;
I_u  = -0.05;
I_v  = -0.3;
tau_e  = 1.0;
tau_i  = 0.5;
tau_0  = 0.2;
eps_c  = 0.01;
T = 0;


par0=cellfun(@(x)evalin('caller',x),parnames(1:ind.tau_0));
for k = 1:N-1
    par0(ind.(sprintf('tau_%d', k))) = 0;
end

ind_taus = arrayfun(@(k) ind.(sprintf('tau_%d',k)), 1:N-1);

%% set user-defined functions
fsymbolic=set_symfuncs(@sym_dWC_sync,'sys_tau',@() [ind.tau_0, ind_taus],'sys_cond',@(pt)cond_dWC_sync(pt,ind));
funcs=fsymbolic;

%% Solve for eps_c=0 using dde23
tspan = [0,300];
history = @(t) [0.4; 0.6];

ddefun = @(t,y,Z) rhs_eps0(t,y,Z,par0,ind);
lags = par0(ind.tau_0);

% Detect upward crossings of u(t) = u_star.
u_star = 0.5;
eventfun = @(t,y,Z) poincare_event(t,y,Z,u_star);
opts = ddeset( ...
    'RelTol', 1e-6, ...
    'AbsTol', 1e-6, ...
    'MaxStep', 0.1, ...
    'Events', eventfun);

% solve
sol = dde23(ddefun, lags, history, tspan, opts);

% Plot solution
tt = linspace(tspan(1), tspan(2), 20000);
yy = deval(sol, tt);

figure(1); clf;
plot(tt, yy(1,:), 'LineWidth', 1.2);
hold on;
plot(tt, yy(2,:), 'LineWidth', 1.2);
yline(u_star, '--');
plot(sol.xe, sol.ye, 'o');
grid on;
xlabel('t');
legend('u(t)', 'v(t)', 'u_*');
title('Delayed Wilson-Cowan system, eps_c = 0');

figure(2); clf;
plot(yy(1,:), yy(2,:), 'LineWidth', 1.2);
hold on;
grid on;
xlabel('u');
ylabel('v');
title('Phase portrait');


%% Setup branch

% If tau_1 is fixed and tau_2,...,tau_10 are constrained by tau_1:
corpar = [ind_taus(2:end), ind.T];

free_pars = [ind_taus, ind.T];

rot_br = branch_from_sol(funcs, sol, free_pars, par0, ...
    'indperiod', ind.T, ...
    'corpar', corpar, ...
    'extra_condition', true, ...
    'print_residual_info', 1);

%% Continue in tau_1 and compute stability

% Bounds and step size for continuation in the base connection delay tau_1.
rot_br.parameter.min_bound = [ind.tau_1, 0.0];
rot_br.parameter.max_bound = [ind.tau_1, 50.0];
rot_br.parameter.max_step  = [ind.tau_1, 0.1];

% Continue forward.
figure(1); clf; hold on;
rot_br.method.continuation.plot = 1;
[rot_br, suc] = br_contn(funcs, rot_br, 1000);
if suc == 0
    warning('Forward continuation stopped before completing all requested steps.');
end
hold off;

% Compute Floquet stability of periodic orbits.
rot_br = br_stabl(funcs, rot_br, 0, 1);
nunst = GetStability(rot_br, 'exclude_trivial', true);

%%
figure(2); clf; hold on;
plot(rot_br.point(1).mesh,rot_br.point(1).profile,'-');
plot(rot_br.point(end).mesh,rot_br.point(end).profile, ':');
hold off;

%% plot bifurcation diagram

tau1_vals = arrayfun(@(p) p.parameter(ind.tau_1), rot_br.point);
per_vals = arrayfun(@(p) p.parameter(ind.T), rot_br.point);

stable = nunst == 0;
unstable = ~stable;

figure(3); clf; hold on;

plot(tau1_vals, per_vals, '-k', ...
    'HandleVisibility', 'off');

plot(tau1_vals(stable), per_vals(stable), '.b', ...
    'MarkerSize', 14, ...
    'DisplayName', 'stable');

plot(tau1_vals(unstable), per_vals(unstable), 'xr', ...
    'MarkerSize', 7, ...
    'DisplayName', 'unstable');

grid on;
xlabel('\tau_1');
ylabel('period');
title('Synchronous branch continued in \tau_1');
legend('Location', 'northwest');

%% save data and figure

% Save MATLAB workspace
save('dWC_tau1_continuation_workspace.mat');

% Save bifurcation diagram as EPS
figure(3);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, 'dWC_tau1_bifurcation_diagram.eps', '-depsc', '-painters');