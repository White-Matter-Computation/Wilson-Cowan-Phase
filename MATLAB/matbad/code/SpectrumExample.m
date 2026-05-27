% Spectrum of linearised problem 

clc; clear all; close all; 
% Loading the path where the DDE-Biftool scripts are stored
addpath('../DDE-Biftool/ddebiftool',...
    '../DDE-Biftool/ddebiftool_extra_psol',...
    '../DDE-Biftool/ddebiftool_extra_nmfm',...
    '../DDE-Biftool/ddebiftool_extra_symbolic',...
    '../DDE-Biftool/ddebiftool_utilities');

% x'(t) = a00 x(t) + a0 x(t-tau0) + a1 x(t-tau1{1}) 
% Avec = [a00, a0, a1, ...];
% Tvec = [tau0, tau1, ...]
%Avec = [1, 2 , 3];
Tau = [0, 1, 3];
%A = reshape(Avec, 1, 1, numel(Avec)); % dim x dim x (m+1)
A(:,:,1) = [[-10,-2];[-3,-40]];
A(:,:,2) = [[5,2];[1,3]];
A(:,:,3) = [[3,8];[6,4]];

lam = dde_stst_eig_cheb(A, Tau, ...
                'minimal_real_part', -2, ...
                'root_accuracy', 1e-8, ...
                'max_number_of_eigenvalues', 20, ...
                'min_number_of_eigenvalues', 10);

%% Plotting
figure(1); clf;
plot(lam.l0, 'oy');