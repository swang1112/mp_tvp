clc
clear all
addpath('fun')

%% model specification
K = 5; % number of factors
T = 100; % number of effective observations

%% Block II: Draw a1, phi and sigma_eta (coproduct: eps_1, omic)


% priors
M_a1 % Kx1
V_a1 % KxK

M_phi % 1x1
V_phi % 1x1

t_eta % 1x1
v_eta % 1x1

% input: 
% Z: Tx1, HFI
% U: TxK, VAR residuals

% 1. Draw a1 | z, y, phi, sigma_eta

omic = phi/sigma_eta; % s2n-ratio based on last draw

V_a1_star = inv(omic^2 * U'*U + inv(V_a1));
M_a1_star = V_a1_star * (omic^2 * U'*Z + inv(V_a1)*M_a1);
a1        = M_a1_star+(randn(1,K)*chol(V_a1_star))';

eps_1 = U*a1; % mp shocks

% 2. Draw phi | z, y, a1, sigma_eta
V_phi_star = 1/(eps_1'*eps_1/sigma_eta^2 + 1/V_phi);
M_phi_star = V_phi_star * (eps_1'*Z/sigma_eta^2 + M_phi/V_phi);
phi        = M_phi_star+(randn(1,K)*chol(V_phi_star))';

eta = Z - phi * eps_1; % measurement error in HFI

% 3. Draw sigma_eta^2 | z, y, a1, phi
t_eta_star = t_eta + T;
v_eta_star = v_eta + eta'*eta;
foo0       = randn(t_eta_star,1);
sigma_eta2 = v_eta_star / foo0'*foo0;
sigma_eta  = sqrt(sigma_eta2);

omic_star  = phi_star / sigma_eta;



