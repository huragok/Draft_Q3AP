clear all;
close all;
clc;

addpath('./functions');

%% 1. Simulation settings
% 32QAM
Eb2N0_noncore = [-3.5, -3, -2.5, -2, -1.5, -1, -0.8, -0.6, -0.4, -0.3, -0.2]; % 1/sigma2 in dB
Eb2N0_seddik = [-4.5, -4, -3.5, -3, -2.6, -2.2, -2, -1.9, -1.8, -1.7, -1.6];
Eb2N0_Q3AP = [-5, -4.5, -4, -3.5, -3.2, -3, -2.9, -2.8, -2.7, -2.6, -2.5];

Nbps = 5;
type_mod = 'QAM';
pwr = 1;
X = get_constellation(Nbps, type_mod, pwr);
Q = 2 ^ Nbps;
K = 10;
amp = 1;

max_frame = 2000;
iter_max = 5;
coding_rate = 3 / 4;
nldpc = 2400;

seed = 8;

map_noncore = [1 : Q; 1 : Q];
map_seddik = [11, 8, 9, 10, 23, 20, 21, 22, 3, 0, 1, 2, 15, 12, 13, 14, 27, 24, 25, 26, 7, 4, 5, 6, 19, 16, 17, 18, 31, 28, 29, 30;
              11, 8, 9, 10, 23, 20, 21, 22, 3, 0, 1, 2, 15, 12, 13, 14, 27, 24, 25, 26, 7, 4, 5, 6, 19, 16, 17, 18, 31, 28, 29, 30] + 1; % A very simple remap based on "Trans-Modulation in Wireless Relay Networks"
map_Q3AP = [30, 10, 12, 15, 8, 14, 31, 19, 23, 20, 22, 2, 25, 16, 24, 5, 13, 6, 29, 27, 21, 26, 28, 7, 0, 4, 1, 18, 17, 11, 9, 3;
            30, 26, 12, 14, 13, 11, 27, 3, 21, 20, 22, 2, 29, 16, 25, 7, 24, 6, 23, 10, 9, 31, 28, 5, 0, 4, 1, 18, 17, 15, 8, 19] + 1;

%% 2. Initialization
n_Eb2N0_noncore = length(Eb2N0_noncore);
n_Eb2N0_seddik = length(Eb2N0_seddik);
n_Eb2N0_Q3AP = length(Eb2N0_Q3AP);

%% 2. Initialization
mu_h = sqrt(K / (K + 1)) * [1; 1; amp];
sigma2_h = 1 / (K + 1) * [1; 1; abs(amp) ^ 2];

%% 3. Not let us test the coded bit error rate
codedBER_noncore = zeros(1, n_Eb2N0_noncore); % coded BER for noncore scheme
codedBER_seddik = zeros(1, n_Eb2N0_seddik);
codedBER_Q3AP = zeros(1, n_Eb2N0_Q3AP);

% The non-CoRe waterfall curve
sigma2_v = 1 ./ (Nbps * 10 .^ (Eb2N0_noncore / 10));
for i_Eb2N0 = 1 : n_Eb2N0_noncore
    tic
    % Compute the bit error rate using Monte-Carlo simulation
    codedBER_noncore(i_Eb2N0) = get_codedBER(X, map_noncore, mu_h, sigma2_h, sigma2_v(i_Eb2N0), max_frame, iter_max, coding_rate, nldpc, seed);
    toc;
    disp(['Non-CoRe, Eb/N0 = ', num2str(Eb2N0_noncore(i_Eb2N0)), 'dB, coded BER = ', num2str(codedBER_noncore(i_Eb2N0))]);
end

% The Seddik waterfall curve
sigma2_v = 1 ./ (Nbps * 10 .^ (Eb2N0_seddik / 10));
for i_Eb2N0 = 1 : n_Eb2N0_seddik
    tic
    % Compute the bit error rate using Monte-Carlo simulation
    codedBER_seddik(i_Eb2N0) = get_codedBER(X, map_seddik, mu_h, sigma2_h, sigma2_v(i_Eb2N0), max_frame, iter_max, coding_rate, nldpc, seed);
    toc;
    disp(['Seddik, Eb/N0 = ', num2str(Eb2N0_seddik(i_Eb2N0)), 'dB, coded BER = ', num2str(codedBER_seddik(i_Eb2N0))]);
end

% The Q3AP waterfall curve
sigma2_v = 1 ./ (Nbps * 10 .^ (Eb2N0_Q3AP / 10));
for i_Eb2N0 = 1 : n_Eb2N0_Q3AP
    tic
    % Compute the bit error rate using Monte-Carlo simulation
    codedBER_Q3AP(i_Eb2N0) = get_codedBER(X, map_Q3AP, mu_h, sigma2_h, sigma2_v(i_Eb2N0), max_frame, iter_max, coding_rate, nldpc, seed);
    toc;
    disp(['Q3AP, Eb/N0 = ', num2str(Eb2N0_Q3AP(i_Eb2N0)), 'dB, coded BER = ', num2str(codedBER_Q3AP(i_Eb2N0))]);
end

%% 4. Visualization
% The coded BER
h = figure;
semilogy(Eb2N0_noncore, codedBER_noncore, 'k+-', 'linewidth', 2), hold on;
semilogy(Eb2N0_seddik, codedBER_seddik, 'b^--', 'linewidth', 2);
semilogy(Eb2N0_Q3AP, codedBER_Q3AP, 'ro-.', 'linewidth', 2);

grid on;
set(gca, 'Fontsize', 18);
xlabel('E_b/N_0(dB)'), ylabel('Coded BER');
legend({'Gray', 'Seddik', 'MoDiv'}, 'Location', 'northeast');
saveas(h, 'waterfall_32QAM.fig');

save(['waterfall_32QAM.mat'], 'Eb2N0_noncore', 'Eb2N0_seddik', 'Eb2N0_Q3AP', 'codedBER_noncore', 'codedBER_seddik', 'codedBER_Q3AP')
