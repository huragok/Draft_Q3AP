clear all;
close all;
clc;

addpath('./functions/')

%% 1. Simulation settings
Nbps = 4;
type_mod = 'QAM';
pwr = 1;
X = get_constellation(Nbps, type_mod, pwr);
Q = 2 ^ Nbps;
K = 10;
amp = 1;
Eb2N0 = [-2 : 4]; % Eb/N0 in dB
n_Eb2N0 = length(Eb2N0);

max_frame = 100;
iter_max = 5;
coding_rate = 3 / 4;
nldpc = 2400;

seed = 8;

map_hans = [13 1 5 9 12 0 4 8 15 3 7 11 14 2 6 10;
            7 3 15 11 4 0 12 8 5 1 13 9 6 2 14 10] + 1;
map_karim = [5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0;
             5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0] + 1; % A very simple remap based on "Trans-Modulation in Wireless Relay Networks"
map_uniform = [1 : Q; 1 : Q];

%% 2. Initialization
mu_h = sqrt(K / (K + 1)) * [1; 1; amp];
sigma2_h = 1 / (K + 1) * [1; 1; abs(amp) ^ 2];
sigma2_v = 1 ./ (Nbps * 10 .^ (Eb2N0 / 10));

%% 3. Compute the BER
BER_coded_uniform = zeros(n_Eb2N0, 1);
BER_coded_hans = zeros(n_Eb2N0, 1);
BER_coded_karim = zeros(n_Eb2N0, 1);

%matlabpool open 4
for i_Eb2N0 = 1 : n_Eb2N0
    
    tic;
    BER_coded_uniform(i_Eb2N0) = get_codedBER(X, map_uniform, mu_h, sigma2_h, sigma2_v(i_Eb2N0), max_frame, iter_max, coding_rate, nldpc, seed);    
    BER_coded_hans(i_Eb2N0) = get_codedBER(X, map_hans, mu_h, sigma2_h, sigma2_v(i_Eb2N0), max_frame, iter_max, coding_rate, nldpc, seed);
    BER_coded_karim(i_Eb2N0) = get_codedBER(X, map_karim, mu_h, sigma2_h, sigma2_v(i_Eb2N0), max_frame, iter_max, coding_rate, nldpc, seed);
   

    disp(['Test case ', num2str(i_Eb2N0)]);
    disp(['Uniform mapping, coded BER = ', num2str(BER_coded_uniform(i_Eb2N0))]);
    disp(['Hans remapping, coded BER = ', num2str(BER_coded_hans(i_Eb2N0))]);
    disp(['Karim remapping, coded BER = ', num2str(BER_coded_karim(i_Eb2N0))]);
    
    toc;
    
end
%matlabpool close


%% 4. Visualization
figure;
semilogy(Eb2N0, BER_coded_uniform, 'bo--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_coded_karim(range), 'rs--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_coded_hans(range), 'm^--', 'linewidth', 2), hold on;
grid on;
set(gca, 'Fontsize', 16);
xlabel('E_b/N_0(dB)'), ylabel('Coded BER'), legend('Uniform', 'Karim', 'Hans');
