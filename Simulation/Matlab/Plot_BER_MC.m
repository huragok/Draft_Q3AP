clear all;
close all;
clc;

load('Rician_16QAM_K10_amp.mat');
BER_uniform_10 = BER_uniform;
BER_seddik_10 = BER_karim;
BER_Q3AP_10 = BER_hans;

load('Rician_16QAM_K5_amp.mat');
BER_uniform_5 = BER_uniform;
BER_seddik_5 = BER_karim;
BER_Q3AP_5 = BER_hans;

amp = [1, sqrt(2), 2];
Eb2N0 = [-2 : 4];

figure;
semilogy(Eb2N0, BER_uniform_10(1 : 7), 'bo:', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_seddik_10(1 : 7), 'bs--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_Q3AP_10(1 : 7), 'b^-', 'linewidth', 2), hold on;

semilogy(Eb2N0, BER_uniform_10(8 : 14), 'ro:', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_seddik_10(8 : 14), 'rs--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_Q3AP_10(8 : 14), 'r^-', 'linewidth', 2), hold on;

semilogy(Eb2N0, BER_uniform_10(15 : 21), 'mo:', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_seddik_10(15 : 21), 'ms--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_Q3AP_10(15 : 21), 'm^-', 'linewidth', 2), hold on;
grid on;
set(gca, 'Fontsize', 18);
xlabel('E_b/N_0(dB)'), ylabel('BER'), ylim([1e-6, 1e-1])
h_legend = legend('Gray, a = 1', 'Seddik, a = 1', 'Q3AP, a = 1',...
    'Gray, a = 2^{1/2}', 'Seddik, a = 2^{1/2}', 'Q3AP, a = 2^{1/2}',...
    'Gray, a = 2', 'Seddik, a = 2', 'Q3AP, a = 2', 'Location','SouthWest');
set(h_legend,'FontSize',16);

figure;
semilogy(Eb2N0, BER_uniform_5(1 : 7), 'bo:', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_seddik_5(1 : 7), 'bs--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_Q3AP_5(1 : 7), 'b^-', 'linewidth', 2), hold on;

semilogy(Eb2N0, BER_uniform_5(8 : 14), 'ro:', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_seddik_5(8 : 14), 'rs--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_Q3AP_5(8 : 14), 'r^-', 'linewidth', 2), hold on;

semilogy(Eb2N0, BER_uniform_5(15 : 21), 'mo:', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_seddik_5(15 : 21), 'ms--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_Q3AP_5(15 : 21), 'm^-', 'linewidth', 2), hold on;
grid on;
set(gca, 'Fontsize', 18);
xlabel('E_b/N_0(dB)'), ylabel('BER'), ylim([1e-6, 1e-1])
h_legend = legend('Gray, a = 1', 'Seddik, a = 1', 'Q3AP, a = 1',...
    'Gray, a = 2^{1/2}', 'Seddik, a = 2^{1/2}', 'Q3AP, a = 2^{1/2}',...
    'Gray, a = 2', 'Seddik, a = 2', 'Q3AP, a = 2', 'Location','SouthWest');
set(h_legend,'FontSize',16);
