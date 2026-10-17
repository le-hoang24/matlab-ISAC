function [D_r, SNRr] = echo_generate(tpara, para, W, X)
% Generate noise-free echo signals for MIMO-OFDM ISAC sensing.
%
% Inputs:
%   tpara : Target parameters struct (theta_rad, R, v, N)
%   para  : System parameters struct (Nr, Ns, Nt, L, dt, dr, fc, c, delta_f, T_s_cp, sigma_s, K)
%   W     : Precoding matrix (Nt x K x Ns)
%   X     : Baseband transmit signal (Nt x L x Ns)
%
% Outputs:
%   D_r  : Noise-free echo signal (Nr x Ns x L)
%   SNRr : Receive SNR (scalar, meaningful only for single-target case)
%
% Note: Angle and range digital frequencies are negated for convenient
%       index mapping in subsequent estimation. Output excludes path loss
%       and RCS effects; noise is added externally when sweeping SNR.

snr_num = 0;

%% Received echo from targets after DFT (frequency domain)
D_r = zeros(para.Nr, para.Ns, para.L);
phase_shift_D = exp(1i*2*pi*para.T_s_cp * (2*tpara.v*para.fc/para.c) * (0:1:para.L-1));

for i = 1:para.Ns
    B_steering = exp(1i*2*pi * ((0:1:para.Nt-1).') * para.dt * sin(tpara.theta_rad.') * (para.fc) / para.c);
    A_steering = exp(1i*2*pi * ((0:1:para.Nr-1).') * para.dr * sin(tpara.theta_rad.') * (para.fc) / para.c);
    Rl = tpara.R + 0*tpara.v * [0:1:para.L-1] * para.T_s_cp;
    phase_shift_R = exp(1i*2*pi * (i-1) * para.delta_f * 2 * (Rl) / para.c);
    temp = reshape(A_steering * (B_steering'*X(:,:,i) .* (phase_shift_R .* phase_shift_D)), [para.Nr, 1, para.L]);
    D_r(:,i,:) = temp;
    if tpara.N == 1
        snr_num = snr_num + B_steering'*W(:,:,i) * W(:,:,i)'*B_steering;
    end
end

%% Receive SNR
if tpara.N == 1
    SNRr = snr_num / (para.Ns * para.sigma_s^2);
else
    SNRr = 1;
end

end
