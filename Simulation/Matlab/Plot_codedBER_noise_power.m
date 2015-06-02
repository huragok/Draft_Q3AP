clear all;
close all;
clc;

addpath('./functions');

%% 1. Simulation settings
% 16QAM
Eb2N0_noncore = [-5, -4.5, -3.5, -3, -2.8, -2.6, -2.4, -2.3, -2.2, -2, -1.9, -1.85]; % 1/sigma2 in dB
Eb2N0_seddik = [-5.5, -5, -4.5, -4.2, -4, -3.8, -3.6, -3.5, -3.4, -3.2, -3.15];
Eb2N0_Q3AP = [-6, -5.5, -5, -4.6, -4.4, -4.2, -4.1, -4, -3.9, -3.8, -3.75];

Nbps = 4;
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
map_seddik = [5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0;
             5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0] + 1; % A very simple remap based on "Trans-Modulation in Wireless Relay Networks"
map_Q3AP = [13 1 5 9 12 0 4 8 15 3 7 11 14 2 6 10;
            7 3 15 11 4 0 12 8 5 1 13 9 6 2 14 10] + 1;

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
saveas(h, 'waterfall_16QAM.fig');

save(['waterfall_16QAM.mat'], 'Eb2N0_noncore', 'Eb2N0_seddik', 'Eb2N0_Q3AP', 'codedBER_noncore', 'codedBER_seddik', 'codedBER_Q3AP')
