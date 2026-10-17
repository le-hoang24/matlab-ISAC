% Fig. 3 baseline: user != target; author's Proposed estimator only.
% Based on Rang Liu's main_fig3.m (MIT, see LICENSE).
% [CHANGE] comments identify modifications relative to the original.
clear; clc; close all;

%% [CHANGE 1] Reproducible, shorter Monte Carlo experiment
seed = 20261009;
rng(seed, 'twister');
root_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(root_dir, 'function'));
Ntest = 20;                          % OLD: 500; quick smoke check: 3
N_trial = 2;                         % OLD: 10; quick smoke check: 1
SNR_dB_range = -60:5:40;             % Unchanged
range_Nr = [4, 16, 64];              % Unchanged (Nt remains 16)
compare_old_geometry = true;        % false: run separated baseline only
angle_offset_deg = 10;              % Separation of communication user/target
% This reproduces the original received-SNR sweep, not a transmit-power sweep.
scenario_names = {'separated'};
if compare_old_geometry
    scenario_names = {'matched', 'separated'};
end
nScenario = numel(scenario_names);
nNr = numel(range_Nr);
nSNR = numel(SNR_dB_range);

%% System parameters: unchanged from original main_fig3.m
para.mod_order = 16;
para.Nt = 16;
para.Ns = 512;
para.L = 256;
para.c = 3e8;
para.fc = 28e9;
para.dt = para.c / para.fc / 2;
para.dr = para.c / para.fc / 2;
para.delta_f = 120e3;
para.T_s = 1 / para.delta_f;
para.T_cp = 0.59e-6;
para.T_s_cp = para.T_s + para.T_cp;
para.K = 1;
tpara.N = 1;                        % Explicitly one target, independently of K
para.sigma_c = sqrt(10^(-90/10));
para.sigma_s = sqrt(10^(-90/10));
a_zero_padding = 1;
tpara.N_fft_d = a_zero_padding * para.Ns;
tpara.N_fft_v = a_zero_padding * para.L;

%% [CHANGE 2] One row per Nr (only Proposed), with separate geometry dimension
RMSE_sum = zeros(nNr, 3, nSNR, nScenario);
target_log = zeros(Ntest, 3);
user_log = zeros(Ntest, 2, nScenario); % [angle(deg), distance(m)]
target_gain = zeros(Ntest, nScenario);

for i_HChange = 1:Ntest
    %% Target generation: unchanged distributions
    tpara.theta_deg = -30 + 60 * rand(tpara.N, 1);
    tpara.theta_rad = deg2rad(tpara.theta_deg);
    tpara.R = 40 + 40 * rand(tpara.N, 1);
    tpara.v = -50 + 100 * rand(tpara.N, 1);
    true_adv = [tpara.theta_deg, tpara.R, tpara.v];
    target_log(i_HChange,:) = true_adv;

    %% [CHANGE 3] User geometry independent of target identity
    % Fixed 10-deg separation, reflected inward at +30 deg to retain the FoV.
    separated_angle_deg = tpara.theta_deg + angle_offset_deg;
    if separated_angle_deg > 30
        separated_angle_deg = tpara.theta_deg - angle_offset_deg;
    end
    % Identical QAM samples in matched/separated cases and for all Nr.
    DATA = randi([0 para.mod_order-1], para.K, para.L, para.Ns);
    S = qammod(DATA, para.mod_order, 'UnitAveragePower', true);
    noise_seed = randi(2^31-1);
    next_geometry_rng = rng;        % Geometry does not depend on trial counts

    for i_case = 1:nScenario
        if strcmp(scenario_names{i_case}, 'matched')
            user_a = tpara.theta_rad;    % OLD: user == target
            user_d = tpara.R;
        else
            user_a = deg2rad(separated_angle_deg);
            user_d = tpara.R + 30;       % Distinct location; not used by LoS H
        end
        user_log(i_HChange,:,i_case) = [rad2deg(user_a), user_d];
        rng(noise_seed, 'twister');      % Paired noise across geometry cases

        %% LoS communication channel: now generated from USER angle
        H = zeros(para.Nt, para.K, para.Ns);
        for i = 1:para.Ns
            for k = 1:para.K
                H(:,k,i) = exp(1i*2*pi * (0:para.Nt-1).' ...
                    * para.dt * sin(user_a(k)) * para.fc / para.c);
            end
        end

        %% W_ZF: formula and Frobenius normalization unchanged
        W = zeros(para.Nt, para.K, para.Ns);
        for i = 1:para.Ns
            W(:,:,i) = H(:,:,i) / (H(:,:,i)' * H(:,:,i));
            W(:,:,i) = 1/norm(W(:,:,i), 'fro') * W(:,:,i);
        end
        X0 = pagemtimes(W, S);
        bt = exp(1i*2*pi*(0:para.Nt-1).' * para.dt ...
            * sin(tpara.theta_rad) * para.fc / para.c);
        target_gain(i_HChange,i_case) = abs(bt' * W(:,:,1))^2;

        for i_Nr = 1:nNr
            para.Nr = range_Nr(i_Nr);
            tpara.N_fft_a = a_zero_padding * para.Nr;

            %% [CHANGE 4] Only author's Proposed; original signal processing
            [D_r0, SNRr0] = echo_generate(tpara, para, W, X0);
            SNRr0 = real(SNRr0);
            assert(isfinite(SNRr0) && SNRr0 > 0, ...
                'Target has zero/invalid illumination under W_ZF.');
            [Y_nfree1_0, sigma_n1] = my_proposed(para, tpara, D_r0, X0);
            assert(all(isfinite(Y_nfree1_0(:))) && all(isfinite(sigma_n1(:))), ...
                'Nonfinite estimator output: inspect near-zero mixed_MS.');

            for i_SNR = 1:nSNR
                fprintf('Geometry %d/%d | Channel %d/%d | Nr=%d | SNR=%g dB\n', ...
                    i_case,nScenario,i_HChange,Ntest,para.Nr,SNR_dB_range(i_SNR));
                SNR = 10^(0.1 * SNR_dB_range(i_SNR));
                para.sigma_beta = sqrt(SNR / SNRr0); % Original normalization
                Y_nfree1 = para.sigma_beta * Y_nfree1_0;

                for i_trial = 1:N_trial
                    % Keep author's processed-noise approximation unchanged.
                    Noise = zeros(tpara.N_fft_a,tpara.N_fft_d,tpara.N_fft_v);
                    for na = 1:tpara.N_fft_a
                        Noise(na,:,:) = sigma_n1(na) .* ...
                            (randn(1,tpara.N_fft_d,tpara.N_fft_v) ...
                            + 1i*randn(1,tpara.N_fft_d,tpara.N_fft_v)) / sqrt(2);
                    end
                    Y = Y_nfree1 + Noise;
                    [estimated_adv, ~] = peak_finding(Y, para, tpara);
                    RMSE_sum(i_Nr,:,i_SNR,i_case) = ...
                        RMSE_sum(i_Nr,:,i_SNR,i_case) + abs(estimated_adv-true_adv).^2;
                end
            end
            clear D_r0 Y_nfree1_0 Y_nfree1 Noise Y;
        end
    end
    rng(next_geometry_rng);           % Resume independent geometry/QAM stream
end
fRMSE = sqrt(RMSE_sum / Ntest / N_trial);

%% [CHANGE 5] Three plots; only Proposed, three Nr curves per geometry
output_dir = fullfile(root_dir, 'baseline_results');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
ylabels = {'RMSE of angle estimation (degree)', ...
    'RMSE of range estimation (m)', 'RMSE of velocity estimation (m/s)'};
metrics = {'angle','range','velocity'};
markers = {'^','s','o'};
colors = lines(nNr);
for i_case = 1:nScenario
    for plot_term = 1:3
        fig = figure('Color','w'); hold on;
        for i_Nr = 1:nNr
            semilogy(SNR_dB_range, reshape(fRMSE(i_Nr,plot_term,:,i_case),1,[]), ...
                ['-' markers{i_Nr}], 'Color',colors(i_Nr,:), 'LineWidth',1.5, ...
                'DisplayName',sprintf('Proposed, Nr = %d',range_Nr(i_Nr)));
        end
        set(gca,'YScale','log'); grid on; legend('Location','best');
        xlabel('Received sensing SNR (dB)'); ylabel(ylabels{plot_term});
        title(sprintf('W ZF / Proposed / %s user-target geometry',scenario_names{i_case}));
        stem = fullfile(output_dir,[scenario_names{i_case} '_' metrics{plot_term}]);
        savefig(fig,[stem '.fig']); exportgraphics(fig,[stem '.png'],'Resolution',160);
    end
end

% Optional overlays make the old/new geometry comparison visible directly.
if compare_old_geometry
    styles = {'--','-'};
    for plot_term = 1:3
        fig = figure('Color','w'); hold on;
        for i_Nr = 1:nNr
            for i_case = 1:nScenario
                semilogy(SNR_dB_range,reshape(fRMSE(i_Nr,plot_term,:,i_case),1,[]), ...
                    [styles{i_case} markers{i_Nr}], 'Color',colors(i_Nr,:), ...
                    'LineWidth',1.5, 'DisplayName',sprintf('%s, Nr = %d', ...
                    scenario_names{i_case},range_Nr(i_Nr)));
            end
        end
        set(gca,'YScale','log'); grid on; legend('Location','best');
        xlabel('Received sensing SNR (dB)'); ylabel(ylabels{plot_term});
        title('Same W ZF rule and Proposed estimator: old/new geometry');
        stem = fullfile(output_dir,['comparison_' metrics{plot_term}]);
        savefig(fig,[stem '.fig']); exportgraphics(fig,[stem '.png'],'Resolution',160);
    end
end
save(fullfile(output_dir,'baseline_results.mat'), 'fRMSE','RMSE_sum', ...
    'SNR_dB_range','range_Nr','scenario_names','Ntest','N_trial','seed', ...
    'angle_offset_deg','target_log','user_log','target_gain','para','tpara');
fprintf('\nDone. Results saved in %s\n',output_dir);
