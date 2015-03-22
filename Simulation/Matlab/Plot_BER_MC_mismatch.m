clear all;
close all;
clc;

Nbps = 4;
type_mod = 'QAM';
pwr = 1;
X = get_constellation(Nbps, type_mod, pwr);
Q = 2 ^ Nbps;
K = 10; % K = 5;
amp = sqrt(2) * exp(1j * pi / 12);
Eb2N0 = [-2 : 4]; % Eb/N0 in dB
N = [100, 200, 200, 300, 400, 400, 400]; % Number of serial expansion
n_amp = length(amp);
n_Eb2N0 = length(Eb2N0);

xi = 1 / 4;
M = 10 ^ 7;

%% 2. Channel settings
channel = 'Rician'; % Channel model, can be specified as AWGN, Rayleigh, Rician, Rayleigh_imp, Rician_imp

map_hans_mismatch = cell(n_Eb2N0, 1);
map_hans_opt = cell(n_Eb2N0, 1);

% Mismatched mapping
map_hans_mismatch{1, 1} = [13 1 5 9 12 0 4 8 15 3 7 11 14 2 6 10;
                  7 3 15 11 4 0 12 8 5 1 13 9 6 2 14 10] + 1;
map_hans_mismatch{2, 1} = [13 6 5 4 3 2 1 0 15 14 7 12 11 10 9 8;
                  13 14 15 12 1 2 3 0 5 6 7 4 9 10 11 8] + 1;
map_hans_mismatch{3, 1} = [15 14 7 12 11 10 9 8 13 6 5 4 3 2 1 0;
                  5 6 7 4 9 10 11 8 13 14 15 12 1 2 3 0] + 1;
map_hans_mismatch{4, 1} = [5 6 7 4 9 10 11 8 13 14 15 12 1 2 3 0;
                  15 14 7 12 11 10 9 8 13 6 5 4 3 2 1 0] + 1;
map_hans_mismatch{5, 1} = [15 1 7 9 12 2 4 10 13 3 5 11 14 0 6 8;
                  15 1 7 9 4 2 12 10 13 3 5 11 6 0 14 8] + 1;
map_hans_mismatch{6, 1} = [5 9 15 3 14 0 4 10 7 11 13 1 12 2 6 8;
                  5 9 15 3 14 0 4 10 7 11 13 1 12 2 6 8] + 1;
map_hans_mismatch{7, 1} = [7 11 13 1 12 2 6 8 5 9 15 3 14 0 4 10;
                  7 11 13 1 12 2 6 8 5 9 15 3 14 0 4 10] + 1;

% The actual mapping
% map_hans_opt{1, 1} = [13 12 15 14 9 8 11 10 5 4 7 6 1 0 3 2;
%                   7 4 5 6 11 8 9 10 15 12 13 14 3 0 1 2] + 1;
% map_hans_opt{2, 1} = [5 9 13 1 6 10 14 2 7 11 15 3 4 8 12 0;
%                   5 11 7 3 14 10 6 2 13 9 15 1 12 8 4 0] + 1;
% map_hans_opt{3, 1} = [15 12 13 14 3 0 1 2 7 4 5 6 11 8 9 10;
%                   15 4 7 6 1 0 3 2 13 12 5 14 9 8 11 10] + 1;
% map_hans_opt{4, 1} = [12 6 5 15 1 8 3 10 13 7 4 14 11 2 9 0;
%                   4 6 5 15 1 8 3 10 13 7 12 14 11 2 9 0] + 1;
% map_hans_opt{5, 1} = [1 6 5 12 11 0 15 10 7 14 3 4 13 8 9 2;
%                   3 6 5 12 11 0 15 10 7 14 1 4 13 8 9 2] + 1;
% map_hans_opt{6, 1} = [9 15 6 3 5 0 4 10 14 11 1 7 12 2 13 8;
%                   9 15 6 3 5 0 4 10 14 11 1 7 12 2 13 8] + 1;
% map_hans_opt{7, 1} = [12 15 1 9 5 10 4 2 3 11 14 13 6 0 7 8;
%                   12 15 1 9 5 10 4 2 3 11 14 13 6 0 7 8] + 1;

              
% The other 2 mappings for comparison              
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
        test_cases(i_case).map_hans_mismatch = map_hans_mismatch{i_Eb2N0, i_amp};
%         test_cases(i_case).map_hans_opt = map_hans_opt{i_Eb2N0, i_amp};
        test_cases(i_case).map_karim  = map_karim;
        test_cases(i_case).map_uniform = map_uniform;
    end
end

n_case = i_case;

%% 3. Compute the BER
BER_uniform = zeros(n_Eb2N0, 1);
BER_hans_mismatch = zeros(n_Eb2N0, 1);
BER_karim = zeros(n_Eb2N0 * n_amp, 1);

BER_uniform_bound = zeros(n_Eb2N0, 1);
BER_hans_mismatch_bound = zeros(n_Eb2N0, 1);
BER_karim_bound = zeros(n_Eb2N0 * n_amp, 1);

for i_case = 1 : n_case
    
    tic;

    BER_uniform(i_case) = get_BER(test_cases(i_case).X, test_cases(i_case).map_uniform, test_cases(i_case).mu_h, test_cases(i_case).sigma2_h, test_cases(i_case).sigma2_eps, test_cases(i_case).sigma2_v, test_cases(i_case).M);
    BER_hans_mismatch(i_case) = get_BER(test_cases(i_case).X, test_cases(i_case).map_hans_mismatch, test_cases(i_case).mu_h, test_cases(i_case).sigma2_h, test_cases(i_case).sigma2_eps, test_cases(i_case).sigma2_v, test_cases(i_case).M);
    BER_karim(i_case) = get_BER(test_cases(i_case).X, test_cases(i_case).map_karim, test_cases(i_case).mu_h, test_cases(i_case).sigma2_h, test_cases(i_case).sigma2_eps, test_cases(i_case).sigma2_v, test_cases(i_case).M);

    BER_uniform_bound(i_case) = get_BER_upper_bound(test_cases(i_case).X, test_cases(i_case).map_uniform, test_cases(i_case).mu_h, test_cases(i_case).sigma2_h, test_cases(i_case).sigma2_eps, test_cases(i_case).sigma2_v, test_cases(i_case).N, test_cases(i_case).xi);
    BER_hans_mismatch_bound(i_case) = get_BER_upper_bound(test_cases(i_case).X, test_cases(i_case).map_hans_mismatch, test_cases(i_case).mu_h, test_cases(i_case).sigma2_h, test_cases(i_case).sigma2_eps, test_cases(i_case).sigma2_v, test_cases(i_case).N, test_cases(i_case).xi);
    BER_karim_bound(i_case) = get_BER_upper_bound(test_cases(i_case).X, test_cases(i_case).map_karim, test_cases(i_case).mu_h, test_cases(i_case).sigma2_h, test_cases(i_case).sigma2_eps, test_cases(i_case).sigma2_v, test_cases(i_case).N, test_cases(i_case).xi);
    
    disp(['Test case ', num2str(i_case)]);

    disp(['Uniform mapping, BER = ', num2str(BER_uniform(i_case))]);
    disp(['Hans mismatched remapping, BER = ', num2str(BER_hans_mismatch(i_case))]);
    disp(['Karim remapping, BER = ', num2str(BER_karim(i_case))]);
    
    disp(['Uniform mapping upperbound, BER = ', num2str(BER_uniform_bound(i_case))]);
    disp(['Hans mismatched remapping upperbound, BER = ', num2str(BER_hans_mismatch_bound(i_case))]);
    disp(['Karim remapping upperbound, BER = ', num2str(BER_karim_bound(i_case))]);
    toc;
    
end

%% 4. Visualization
figure;
semilogy(Eb2N0, BER_uniform_bound, 'bo--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_uniform, 'bo-', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_karim_bound, 'rs--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_karim, 'rs-', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_hans_mismatch_bound, 'm^--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_hans_mismatch, 'm^-', 'linewidth', 2), hold on;
grid on;
set(gca, 'Fontsize', 18);
xlabel('E_b/N_0(dB)'), ylabel('BER');
h_legend = legend('Gray bound', 'Gray MC', 'Seddik bound', 'Seddik MC', 'Q3AP mis bound', 'Q3AP mis MC');
set(h_legend,'FontSize',16);