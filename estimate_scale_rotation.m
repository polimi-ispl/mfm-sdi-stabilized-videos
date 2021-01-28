function [est_scale, est_theta] = estimate_scale_rotation(est_shift, ...
    F_F_noise_logp_deltarho, F_K_logp_deltarho, exp_x_logp_deltarho, F_noise_if_logp_obj)
% [est_scale, est_theta] = estimate_scale_rotation(est_shift, F_F_noise_logp_deltarho, 
% F_logp_K, exp_x_logp_deltarho, Fpolar_obj) estimates the scaling and
% rotation between the noise residual and the shifted PRNU
% INPUTS: 
% est_shift = shift estimated by ga 
% F_F_noise_logp_deltarho = Fourier transform of the MFM-deltarho transform
% of the noise residual
% F_K_logp_deltarho = MFM-deltarho transform of the PRNU
% exp_x_logp_deltarho = exponential term for realigning the PRNU to the
% noise residual
% F_noise_if_logp_obj = output of the function 'LogPolar_tx(F_noise_if)'
% OUTPUTS: 
% est_scale = estimated scaling factor
% est_theta = estimated rotation angle

%% realign noise residual and PRNU according to the estimated shift (in
% Fourier-domain)

% exponential contribution given by the selected shift
exp_contrib_shift = (exp_x_logp_deltarho.^(1i*est_shift(1))).*(circshift(exp_x_logp_deltarho, ...
    [0, floor(size(exp_x_logp_deltarho, 2)/2)]).^(1i*est_shift(2)));
exp_contrib_shift(isnan(exp_contrib_shift)) = 0;

% shift the PRNU
F_F_logp_K_shifted = fft2(F_K_logp_deltarho.*exp_contrib_shift);

%% phase correlation

phi_aux = F_F_noise_logp_deltarho .* conj(F_F_logp_K_shifted);
phi = abs(ifft2(phi_aux ./ abs(eps+phi_aux)));

%% find scale and rotation

est_translation = find_translation_phase_corr(phi); 

thetaIntrinsic = abs(est_translation(1))+1;
logRhoIntrinsic   = abs(est_translation(2))+1;
theta = (thetaIntrinsic-1) .* F_noise_if_logp_obj.thetaMax ./ F_noise_if_logp_obj.numSamplesTheta;
deltaLogRho = (F_noise_if_logp_obj.logRhoMax) ./ F_noise_if_logp_obj.numSamplesRho;
rho = exp((logRhoIntrinsic-1) .* deltaLogRho);

est_theta = -sign(est_translation(1))*theta;
est_scale = rho .^ -sign(est_translation(2));

end

