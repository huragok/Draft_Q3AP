clear all;
close all;
clc;

addpath('./functions/')
%% 1. Generate the Gray mapped constellation
Nbps = 4;
type_mod = 'QAM';
pwr = 1;
X = get_constellation(Nbps, type_mod, pwr);
Q = 2 ^ Nbps;
K = 10;
amp = [1, sqrt(2), 2];
Eb2N0 = [-2 : 4]; % Eb/N0 in dB
N = [100, 200, 200, 300, 400, 400, 400]; % Number of serial expansion
n_amp = length(amp);
n_Eb2N0 = length(Eb2N0);

xi = 1 / 4;
M = 10 ^ 7;

%% 2. Channel settings
channel = 'Rician'; % Channel model, can be specified as AWGN, Rayleigh, Rician, Rayleigh_imp, Rician_imp

map_hans = cell(n_Eb2N0, n_amp);
map_hans{1, 1} = [15 11 7 3 12 8 4 0 13 9 5 1 14 10 6 2;
                  5 9 13 1 4 8 12 0 7 11 15 3 6 10 14 2] + 1;
map_hans{2, 1} = [7 3 15 11 6 2 14 10 5 1 13 9 4 0 12 8;
                  13 1 5 9 14 2 6 10 15 3 7 11 12 0 4 8] + 1;
map_hans{3, 1} = [5 4 13 6 1 0 3 2 7 12 15 14 9 8 11 10;
                  15 12 13 14 3 0 1 2 7 4 5 6 11 8 9 10] + 1;
map_hans{4, 1} = [7 4 5 6 11 8 9 10 15 12 13 14 3 0 1 2;
                  13 12 5 14 9 8 11 10 15 4 7 6 1 0 3 2] + 1;
map_hans{5, 1} = [5 1 7 9 4 0 12 8 13 3 15 11 6 2 14 10;
                  15 3 7 11 12 0 4 8 13 1 5 9 14 2 6 10] + 1;
map_hans{6, 1} = [5 14 7 12 9 2 11 0 15 4 13 6 1 8 3 10;
                  5 14 7 12 11 2 9 0 15 4 13 6 3 8 1 10] + 1;
map_hans{7, 1} = [7 11 13 1 12 2 6 8 5 9 15 3 14 0 4 10;
                  7 11 13 1 12 2 6 8 5 9 15 3 14 0 4 10] + 1;

map_hans{1, 2} = [] + 1;
map_hans{2, 2} = [] + 1;
map_hans{3, 2} = [] + 1;
map_hans{4, 2} = [] + 1;
map_hans{5, 2} = [] + 1;
map_hans{6, 2} = [] + 1;
map_hans{7, 2} = [] + 1;

map_hans{1, 3} = [5 7 11 10 8 9 14 6 13 15 3 2 4 12 0 1;
                  7 5 11 10 9 8 6 14 15 13 3 2 12 4 0 1] + 1;
map_hans{2, 3} = [15 5 12 8 3 1 2 6 7 13 9 11 0 4 10 14;
                  13 15 12 8 1 3 2 6 7 5 9 11 4 0 10 14] + 1;
map_hans{3, 3} = [15 7 9 8 4 12 0 1 5 13 11 10 14 6 2 3;
                  15 7 9 8 12 4 0 1 5 13 11 10 6 14 2 3] + 1;
map_hans{4, 3} = [15 13 14 6 0 1 3 2 5 7 4 12 10 11 9 8;
                  13 5 14 6 0 1 3 2 7 15 4 12 10 11 9 8] + 1;
map_hans{5, 3} = [15 13 3 1 10 14 8 12 5 7 11 9 2 6 0 4;
                  13 5 3 1 10 14 8 12 7 15 11 9 2 6 0 4] + 1;
map_hans{6, 3} = [7 3 13 9 10 14 0 4 15 11 5 1 2 6 8 12;
                  7 3 13 9 10 14 0 4 15 11 5 1 2 6 8 12] + 1;
map_hans{7, 3} = [6 7 1 13 11 10 4 12 15 14 5 0 3 2 9 8;
                  6 7 1 13 11 10 4 12 15 14 5 0 3 2 9 8] + 1;
map_karim = [5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0;
             5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0] + 1; % A very simple remap based on "Trans-Modulation in Wireless Relay Networks"
map_uniform = [1 : Q; 1 : Q];

test_cases = struct();
i_case = 0;
for i_amp = 1 : n_amp
    for i_Eb2N0 = 1 : n_Eb2N0
        i_case = i_case + 1;
        test_cases(i_case).type = 'Rician';
        test_cases(i_case).mu_h = sqrt(K / (K + 1)) * [1; 1; amp(i_amp)];

        test_cases(i_case).sigma2_h = 1 / (K + 1) * [1; 1; abs(amp(i_amp)) ^ 2];
        test_cases(i_case).sigma2_eps = zeros(3, 1);
        test_cases(i_case).Nbps = Nbps;
        test_cases(i_case).X = X;
        test_cases(i_case).Eb2N0 = Eb2N0(i_Eb2N0);
        test_cases(i_case).sigma2_v = 1 / (Nbps * 10 ^ (Eb2N0(i_Eb2N0) / 10));
        test_cases(i_case).N = N(i_Eb2N0);
        test_cases(i_case).xi = xi;
        test_cases(i_case).M = M;
        test_cases(i_case).map_hans = map_hans{i_Eb2N0, i_amp};
        test_cases(i_case).map_karim  = map_karim;
        test_cases(i_case).map_uniform = map_uniform;
    end
end
n_case = i_case;
%% 3. Compute the BER
BER_upperbound_uniform = zeros(n_Eb2N0 * n_amp, 1);
BER_uniform = zeros(n_Eb2N0 * n_amp, 1);
BER_upperbound_hans = zeros(n_Eb2N0 * n_amp, 1);
BER_hans = zeros(n_Eb2N0 * n_amp, 1);
BER_upperbound_karim = zeros(n_Eb2N0 * n_amp, 1);
BER_karim = zeros(n_Eb2N0 * n_amp, 1);
%matlabpool open 4
for i_case = 1 : 7
    
    tic;
    BER_upperbound_uniform(i_case) = get_BER_upper_bound(test_cases(i_case).X, test_cases(i_case).map_uniform, test_cases(i_case).mu_h, test_cases(i_case).sigma2_h, test_cases(i_case).sigma2_eps, test_cases(i_case).sigma2_v, test_cases(i_case).N, test_cases(i_case).xi);
    BER_uniform(i_case) = get_BER(test_cases(i_case).X, test_cases(i_case).map_uniform, test_cases(i_case).mu_h, test_cases(i_case).sigma2_h, test_cases(i_case).sigma2_eps, test_cases(i_case).sigma2_v, test_cases(i_case).M);
    BER_upperbound_hans(i_case) = get_BER_upper_bound(test_cases(i_case).X, test_cases(i_case).map_hans, test_cases(i_case).mu_h, test_cases(i_case).sigma2_h, test_cases(i_case).sigma2_eps, test_cases(i_case).sigma2_v, test_cases(i_case).N, test_cases(i_case).xi);
    BER_hans(i_case) = get_BER(test_cases(i_case).X, test_cases(i_case).map_hans, test_cases(i_case).mu_h, test_cases(i_case).sigma2_h, test_cases(i_case).sigma2_eps, test_cases(i_case).sigma2_v, test_cases(i_case).M);
    BER_upperbound_karim(i_case) = get_BER_upper_bound(test_cases(i_case).X, test_cases(i_case).map_karim, test_cases(i_case).mu_h, test_cases(i_case).sigma2_h, test_cases(i_case).sigma2_eps, test_cases(i_case).sigma2_v, test_cases(i_case).N, test_cases(i_case).xi);
    BER_karim(i_case) = get_BER(test_cases(i_case).X, test_cases(i_case).map_karim, test_cases(i_case).mu_h, test_cases(i_case).sigma2_h, test_cases(i_case).sigma2_eps, test_cases(i_case).sigma2_v, test_cases(i_case).M);

    disp(['Test case ', num2str(i_case)]);
    disp(['Uniform mapping, BER upperbound = ', num2str(BER_upperbound_uniform(i_case))]);
    disp(['Uniform mapping, BER = ', num2str(BER_uniform(i_case))]);
    disp(['Hans remapping, BER upperbound = ', num2str(BER_upperbound_hans(i_case))]);
    disp(['Hans remapping, BER = ', num2str(BER_hans(i_case))]);
    disp(['Karim remapping, BER upperbound = ', num2str(BER_upperbound_karim(i_case))]);
    disp(['Karim remapping, BER = ', num2str(BER_karim(i_case))]);
    toc;
    
    plot_mapping(test_cases(i_case).X, test_cases(i_case).map_hans(1, :), ['Case ', num2str(i_case), 'S-D']);
    plot_mapping(test_cases(i_case).X, test_cases(i_case).map_hans(2, :), ['Case ', num2str(i_case), 'R-D']);
end
%matlabpool close
plot_mapping(X, map_uniform(1, :), 'Gray S-D');
plot_mapping(X, map_uniform(2, :), 'Gray R-D');
save('Rician_16QAM_K10_amp.mat', 'test_cases', 'BER_upperbound_uniform', 'BER_uniform', 'BER_upperbound_hans', 'BER_hans', 'BER_upperbound_karim', 'BER_karim');
%% 4. Visualization
% amp = 1
range = 1 : 7;
figure('Name','Amp = 1')
semilogy(Eb2N0, BER_upperbound_uniform(range), 'bo--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_uniform(range), 'bo-', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_upperbound_karim(range), 'rs--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_karim(range), 'rs-', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_upperbound_hans(range), 'm^--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_hans(range), 'm^-', 'linewidth', 2), hold on;
grid on;
set(gca, 'Fontsize', 16);
xlabel('E_b/N_0(dB)'), ylabel('BER'), legend('Uniform bound', 'Uniform sim.', 'Karim bound', 'Karim sim.', 'Hans bound', 'Hans sim.');

% amp = sqrt(2)
range = 8 : 14;
figure('Name','Amp = sqrt(2)')
semilogy(Eb2N0, BER_upperbound_uniform(range), 'bo--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_uniform(range), 'bo-', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_upperbound_karim(range), 'rs--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_karim(range), 'rs-', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_upperbound_hans(range), 'm^--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_hans(range), 'm^-', 'linewidth', 2), hold on;
grid on;
set(gca, 'Fontsize', 16);
xlabel('E_b/N_0(dB)'), ylabel('BER'), legend('Uniform bound', 'Uniform sim.', 'Karim bound', 'Karim sim.', 'Hans bound', 'Hans sim.');

% amp = 2
range = 15 : 21;
figure('Name','Amp = 2')
semilogy(Eb2N0, BER_upperbound_uniform(range), 'bo--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_uniform(range), 'bo-', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_upperbound_karim(range), 'rs--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_karim(range), 'rs-', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_upperbound_hans(range), 'm^--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_hans(range), 'm^-', 'linewidth', 2), hold on;
grid on;
set(gca, 'Fontsize', 16);
xlabel('E_b/N_0(dB)'), ylabel('BER'), legend('Uniform bound', 'Uniform sim.', 'Karim bound', 'Karim sim.', 'Hans bound', 'Hans sim.');