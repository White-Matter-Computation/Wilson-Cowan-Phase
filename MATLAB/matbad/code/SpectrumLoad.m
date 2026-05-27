% Spectrum of linearised problem 

clc; clear all; close all; 
% Loading the path where the DDE-Biftool scripts are stored
addpath('../DDE-Biftool/ddebiftool',...
    '../DDE-Biftool/ddebiftool_extra_psol',...
    '../DDE-Biftool/ddebiftool_extra_nmfm',...
    '../DDE-Biftool/ddebiftool_extra_symbolic',...
    '../DDE-Biftool/ddebiftool_utilities');

%% Load coefficients from Julia
data  = load('dde_coeffs.mat');
A     = data.A_coeff;   % N x N x (m+1), matches DDEBiftool convention
Tau   = data.Tau;       % [0, tau_0, theta_1, theta_2, ...]

fprintf('System size N = %d\n', size(A,1));
fprintf('Number of delays (incl. 0): %d\n', numel(Tau));
fprintf('Delays: '); disp(Tau);

%% Compute eigenvalues at steady state x=0
lam = dde_stst_eig_cheb(A, Tau, ...
                'minimal_real_part', -2, ...
                'root_accuracy', 1e-8, ...
                'max_number_of_eigenvalues', 20, ...
                'min_number_of_eigenvalues', 10);

%% Plot spectrum
figure(1); clf;
plot(real(lam.l0), imag(lam.l0), 'ob', 'MarkerFaceColor', 'b');
xline(0, '--r'); yline(0, '--k');
xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
title(sprintf('Spectrum (N=%d, %d delays)', size(A,1), numel(Tau)-1));
grid on;