%% Modified Fourier-Mellin code
% a single video frame is tested vs a single PRNU
% Reference:"A MODIFIED FOURIER-MELLIN APPROACH FOR SOURCE DEVICE IDENTIFICATION ON STABILIZED VIDEOS"
% Sara Mandelli, Fabrizio Argenti, Paolo Bestagini, Massimo Iuliani, Alessandro Piva, Stefano Tubaro
% IEEE International Conference on Image Processing (ICIP) 2020
% Requirements: MATLAB Image Processing Toolbox, Global Optimization
% Toolbox, Camera-fingerprint package from "http://dde.binghamton.edu/download/camera_fingerprint"
% (Run the function "compile.m" in folder "CameraFingerprint/Filter").
% @author: Sara Mandelli - sara.mandelli@polimi.it

close all
clearvars
clc

addpath(genpath('CameraFingerprint'));

%%% random seed

rng(42);

%% Genetic algorithm parameters

% number of variables to estimate
nvars = 2;

% options for the optimization
options_ga = optimoptions('ga', 'UseParallel', true, 'MaxStallGenerations',...
    20, 'MaxGenerations', 50, 'PopulationSize', 50);

% lower and upper bounds
tx_min = -90;
tx_max = 90;
ty_min = -90;
ty_max = 90;
lower_b = [tx_min; ty_min];
upper_b = [tx_max; ty_max];

%% other parameters

% frame size for Full-HD sequences
M = 1080;
N = 1920;

% max number of tested frames
max_frames = 10;

% samples for computing the 2D Fourier transform 
F_tx_samples = 2^12;

% angular range for computing the Fourier-Mellin transform
thetaRange = [0 pi];

% exp_x_logp is needed to perform the intial shift correction (see Fig.1 in the original paper) 
x = repmat(1:F_tx_samples, F_tx_samples, 1);
x = (x - 1)./F_tx_samples;
exp_x = fftshift(exp(2*pi*x));
exp_x_logp = LogPolar_tx(exp_x, thetaRange);
exp_x_logp = exp_x_logp.resampledImage; % dimensions: rho x theta

% define number of delta_rho samples to cut the MFM transform (delta_rho = 2 * delta_rho_half)
delta_rho_half = 400; % set it to 400 in order to have delta_rho = 800 
  
%% load device reference fingerprint, scaled and cropped as reported in [1]

load('K.mat');

%% load the i-frame noise residuals from the selected sequence

load('test_noises.mat');
n_frames = size(frame_noises, 3);

%% compute the Modified Fourier Mellin tx for the PRNU

% Fourier-tx of PRNU
F_K = fftshift(fft2(K, F_tx_samples, F_tx_samples));
% convert to log-polar domain
F_K_logp = LogPolar_tx(F_K, thetaRange);
F_K_logp = F_K_logp.resampledImage; % dimensions: rho x theta

% compute the energy of F_K_logp as a function of the rho coordinate
energy_K = sum(abs(F_K_logp).^2, 2);

% cut low frequencies
energy_K = energy_K(1500:end);

% find the rho coordinate corresponding to the peak of the energy
[~, max_idx] = max(energy_K);
crop_rho_right = 1500 + max_idx + delta_rho_half;
crop_rho_left = 1500 + max_idx - delta_rho_half;
if crop_rho_right >= size(F_K_logp, 1)
    crop_rho_right = size(F_K_logp, 1);
end

% crop F_K_logp along the rho coordinate
F_K_logp_deltarho = F_K_logp(crop_rho_left:crop_rho_right, :);

% crop also exp_x_logp in this range of samples along the rho coordinate
exp_x_logp_deltarho = exp_x_logp(crop_rho_left:crop_rho_right, :);

%% select max_frames for the test

% do not consider first frame
test_frames = setdiff(randperm(n_frames, max_frames), 1);

% loop until you select max_frames
while length(test_frames) < max_frames
    
    test_frames =  setdiff(randperm(size(frame_noises, 3), max_frames), 1);
    
end

%% test the selected frames

pce_MFM_deltarho_frames = zeros(max_frames, 1);
pce_onlyshift_frames = zeros(max_frames, 1);

cnt_frames = 1;

for f = test_frames

    noise_frame = frame_noises(:, :, f);
    
    %% compute the PCE without optimizing the process (i.e., search only for a potential shift)
    
    C = crosscorr(noise_frame, K);
    detection = PCE(C, size(noise_frame) - 1);
    pce_onlyshift_frames(cnt_frames) = detection.PCE;
    
    %% run the Modified-Fourier-Mellin approach for testing the frame
    
%     estimate the similarity-transformation parameters (scale,
%     theta, shift_x and shift_y) between the PRNU and the noise residual
%     and compute the PCE between the warped PRNU and the noise residual
    [scale, theta, shift_x, shift_y, pce_frame] = PCE_MFM_deltarho(F_K_logp_deltarho,...
        noise_frame, F_tx_samples, options_ga, lower_b, upper_b, nvars, crop_rho_left, ...
        crop_rho_right, exp_x_logp_deltarho, K);   
    
    % save the PCE value
    pce_MFM_deltarho_frames(cnt_frames) = pce_frame;
    
    cnt_frames = cnt_frames + 1;   
    
    
end

pce_MFM_deltarho_video = max(pce_MFM_deltarho_frames);
fprintf('PCE of the query video = %3.3f\n', pce_MFM_deltarho_video);

figure;
leg = {};
stem(pce_onlyshift_frames);
leg{1} = '$\mathrm{Original \, PCE}$';
hold on,
stem(pce_MFM_deltarho_frames, '*');
leg{2} = '$\mathrm{PCE \, after\,  MFM_{\Delta_{\rho}} \, transform}$';
l = legend(leg);
set(l,'interpreter','latex','fontsize',18,'location','northwest', ...
    'EdgeColor', 0.7*[1, 1, 1]);
grid;
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 12);
xlabel('Frame index','Interpreter','latex','fontsize',18);
xlabel('PCE','Interpreter','latex','fontsize',18);


