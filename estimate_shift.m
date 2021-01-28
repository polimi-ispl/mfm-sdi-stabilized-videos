function phi_out = estimate_shift(shift, F_F_noise_logp_deltarho, F_K_logp_deltarho, exp_x_logp_deltarho)
% phi_out = estimate_shift(shift, F_F_noise_logp_deltarho, F_logp_K, exp_x_logp_deltarho)
% compute the opposite of the maximum of the phase correlation between the
% MFM-deltarho tx of the shifted PRNU and the MFM-deltarho tx of the noise residual
% INPUTS: 
% shift = shift to be applied to the PRNU
% F_F_noise_logp_deltarho = Fourier tx of the MFM-deltarho tx of the noise
% residual
% F_K_logp_deltarho = MFM-deltarho tx of the PRNU K
% exp_x_logp_deltarho = exponential term needed to realign PRNU and noise
% residual
% OUTPUTS:
% phi_out = opposite of the maximum value of the phase correlation

%% realign noise residual and PRNU according to the shift (in Fourier-domain)

% exponential contribution given by the selected shift
exp_contrib_shift = (exp_x_logp_deltarho.^(1i*shift(1))).*(circshift(exp_x_logp_deltarho, ...
    [0, floor(size(exp_x_logp_deltarho, 2)/2)]).^(1i*shift(2)));
exp_contrib_shift(isnan(exp_contrib_shift)) = 0;

% shift the PRNU and take the Fourier tx
F_F_logp_K_shifted = fft2(F_K_logp_deltarho.*exp_contrib_shift);

%% phase correlation

phi_aux = F_F_noise_logp_deltarho .* conj(F_F_logp_K_shifted);
phi = abs(ifft2(phi_aux ./ abs(eps+phi_aux)));

% we can reduce the search space to speed up the process
phi_red = [phi(1:4, 1:70), phi(1:4, end-70:end);
        phi(end-4:end, 1:70), phi(end-4:end, end-70:end)];

% search the maximum of phase correlation
phi_out = -max(phi_red(:));


end
