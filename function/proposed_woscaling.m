function [A, sigma_n] = proposed_woscaling(para, tpara, D_r, X)
% Proposed estimation method without noise scaling compensation.
%
% Inputs:
%   para  : System parameters struct
%   tpara : Target parameters struct
%   D_r   : Echo signal (Nr x Ns x L)
%   X     : Baseband transmit signal (Nt x L x Ns)
%
% Outputs:
%   A       : 3D angle-range-velocity cube after processing
%   sigma_n : Noise scaling factor from processing

%% Step 1: Spatial DFT processing
A = fft(D_r, tpara.N_fft_a, 1);
A = fftshift(A, 1);

%% Step 2: Signal-dependent coefficient removal without scaling
sigma_n = zeros(tpara.N_fft_a, 1);

for i = 1:para.Ns
    a_all = asin((-floor(tpara.N_fft_a/2):1:floor((tpara.N_fft_a-1)/2)).' ...
        * 2*pi/tpara.N_fft_a * para.c / (2*pi*para.dr*(para.fc)));
    B_steering = exp(1i*2*pi * (0:1:para.Nt-1).' * para.dt * sin(a_all.') * (para.fc) / para.c);
    mixed_MS = reshape(B_steering'*X(:,:,i), [tpara.N_fft_a, 1, para.L]);
    A(:,i,:) = A(:,i,:) ./ mixed_MS;
    sigma_n = sigma_n + reshape(sum(abs(1/mixed_MS).^2, [2,3]), tpara.N_fft_a, 1);
end

sigma_n = sqrt(sigma_n * para.Nr * para.sigma_s^2);

%% Step 3: 2D-DFT in fast-time and slow-time dimensions
A = fft(A, tpara.N_fft_d, 2);
A = fft(A, tpara.N_fft_v, 3);
A = fftshift(A, 3);
A = abs(A);

end
