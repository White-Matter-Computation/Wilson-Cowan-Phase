% gen_sym_dWC_sync.m
% Generate symbolic RHS and derivatives for delayed WC system resctricted
% to sync manifold
%
% System:
% tau_e du/dt = -u + f(I_u + w_ee u(t-tau_0)
%                          - w_ei v(t-tau_0)
%                          + eps_c sum_{j=1}^{N-1} w_j u(t-tau_j))
%
% tau_i dv/dt = -v + f(I_v + w_ie u_i(t-tau_0)
%                          - w_ii v_i(t-tau_0))
%
% f(x) = 1/(1+exp(-beta*x))

clear all; close all; clc;

addpath('../../ddebiftool', ...
        '../../ddebiftool_extra_symbolic');

%% State vector

N = 11;         % "network size" 
nvar = 2;       % number of variables
ntau = N;     % number of unique delays [tau0,tau1,...,tauN-1]

x = sym('x', [2, ntau+1]);  % x(:,1) = current state
                            % x(:,2) = t - tau_0
                            % x(:,2+k) = t - tau_k, k = 1,...,floor(N/2)
                                           
%% Parameters
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
    'T', ...
    'tau_0'};

for k = 1:N-1
    parnames{end+1} = sprintf('tau_%d', k); 
end

cind = [parnames; num2cell(1:length(parnames))];
ind = struct(cind{:}); 

syms(parnames{:});
par = sym(parnames);

%% Nonlinearity and weights fixed
f = @(z,beta)  1 ./ ( 1 + exp(-20*z)); % beta = 20;
w = @(k,N)  exp(-min(abs(k), N - abs(k)) / 2); % a = 2;

%% RHS
F = sym(zeros(nvar,1));

coupling = sym(0);
for k = 1:N-1
        wk = w(k,N);
        coupling = coupling + wk*x(1,2+k);
end

z_u = I_u ...
    + w_ee*x(1,2) ...
    - w_ei*x(2,2) ...
    + eps_c*coupling;
z_v = I_v ...
    + w_ie*x(1,2) ...
    - w_ii*x(2,2);

F(1) = (-x(1,1) + f(z_u)) / tau_e;  % u component
F(2) = (-x(2,1) + f(z_v)) / tau_i;  % v component

%% Generate DDE-Biftool functions
[fstr, derivs] = dde_sym2funcs( ...
    F, x, par, ...
    'filename', 'sym_dWC_sync', ...
    'directional_derivative', true);