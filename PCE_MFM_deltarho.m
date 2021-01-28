function [est_scale, est_theta, est_shift_x, est_shift_y, pce_frame] = PCE_MFM_deltarho(F_K_logp_deltarho,...
    noise_frame, F_tx_samples, options_ga, lower_b, upper_b, nvars, ...
    crop_rho_left, crop_rho_right, exp_x_logp_deltarho, K)     
% [est_scale, est_theta, est_shift_x, est_shift_y, pce_frame] = PCE_MFM_deltarho(
% F_K_logp_deltarho, noise_frame, F_tx_samples, options_ga, lower_b, upper_b, 
% nvars, crop_rho_left, crop_rho_right, exp_x_logp_deltarho, K) estimates the 
% similarity transformation between the PRNU K and the noise residual of the selected frame,
% and computes the PCE between the warped PRNU and the noise residual.
% INPUTS: 
% F_K_logp_deltarho = MFM-deltarho tx of the PRNU K
% noise_frame = noise residual of the selected frame
% F_tx_samples = number of Fourier tx samples
% options_ga = options of ga optimizer
% lower_b = lower bound for searching the shift
% upper_b = upper bound for searching the shift
% nvars = number of variables to estimate with ga
% crop_rho_left = left-end point to crop the MFM tx
% crop_rho_right = right-end point to crop the MFM tx
% exp_x_logp_deltarho = exponential term needed to realign PRNU and noise
% residual
% K = the PRNU
% OUTPUTS:
% est_scale = estimated scaling factor
% est_theta = estimated rotation angle
% est_shift_x = estimated shift along x coordinate
% est_shift_y = estimated shift along y coordinate
% pce_frame = PCE between the warped PRNU and the noise residual

%% MFM-deltarho tx of the noise residual

F_noise_if = fftshift((fft2(noise_frame, F_tx_samples, F_tx_samples)));
thetaRange = [0 pi];
F_noise_if_logp_obj = LogPolar_tx(F_noise_if, thetaRange);
F_noise_logp = F_noise_if_logp_obj.resampledImage;
F_noise_logp_deltarho = F_noise_logp(crop_rho_left:crop_rho_right, :);

%% Estimate the shift between the noise residual and the PRNU (Genetic Algorithm)

% compute the Fourier-tx of F_noise_logp_deltarho
F_F_noise_logp_deltarho = fft2(F_noise_logp_deltarho);

% estimate the shift with the genetic algorithm
func = @(x) estimate_shift(x, F_F_noise_logp_deltarho, F_K_logp_deltarho, exp_x_logp_deltarho);
[est_shift, ~] = ga(func, nvars,[], [], [], [], lower_b, upper_b, [], [1 2], options_ga);

%% Realign the PRNU to the noise residual according to the estimated shift

K_realigned = imtranslate(K, [-est_shift(1), est_shift(2)]);

%% Estimate the scale and rotation with phase correlation

[est_scale, est_theta] = estimate_scale_rotation(est_shift, ...
    F_F_noise_logp_deltarho, F_K_logp_deltarho, exp_x_logp_deltarho, F_noise_if_logp_obj);

%% Correct the final estimation of the shift

[est_shift_x, est_shift_y] = solve_for_translation_given_scale_rotation(K, noise_frame, est_scale, est_theta, F_tx_samples);

%% Warp the PRNU according to the estimated scale and rotation

ref = imref2d(size(noise_frame));

t1 = [cos(est_theta), -sin(est_theta), 0;...
    sin(est_theta), cos(est_theta), 0;...
    0, 0, 1];
t2 = [est_scale, 0, 0;...
    0, est_scale, 0; ...
    0, 0, 1];
transf = t2*t1;
tform = affine2d(transf);

[K_warped, ~] = imwarp(K_realigned, tform,'OutputView', ref, 'interp', 'cubic');

%% Compute the PCE between the warped PRNU and the noise residual

C = crosscorr(noise_frame, K_warped);
detection = PCE(C, size(noise_frame) - 1);
pce_frame = detection.PCE;

end