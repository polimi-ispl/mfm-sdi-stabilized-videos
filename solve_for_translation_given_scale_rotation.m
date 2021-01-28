function [est_shift_x, est_shift_y] = solve_for_translation_given_scale_rotation(K,...
    noise_frame, scale, theta, F_tx_samples)
% Inspired to the MATLAB function 'solveForTranslationGivenScaleAndRotation.m'
% There is a 180 degree ambiguity in theta solved in R,Theta space. This
% ambiguity stems from the conjugate symmetry of the Fourier spectrum for real
% valued input images.
%
% This function resolves the ambiguity by forming two resampled versions of K
% rotated by theta, theta+180, phase correlating each version of the
% resampled image with noise_frame, and choose the scale,Theta that has the highest
% final peak correlation during recovery of translation.
% INPUTS: 
% K = device PRNU
% noise_frame = noise residual
% scale = estimated scale as output of the function
% 'estimate_scale_rotation.m'
% theta = estimated rotation as output of the function
% 'estimate_scale_rotation.m'
% F_tx_samples = number of Fourier-tx samples
% OUTPUTS: 
% est_shift_x = final estimated shift along x coordinate
% est_shift_y = final estimated shift along y coordinate

theta1 = theta;
theta2 = theta+pi;

tform1 = affine2d([scale.*cos(theta1) -scale.*sin(theta1) 0; scale.*sin(theta1) scale.*cos(theta1) 0; 0 0 1]);
tform2 = affine2d([scale.*cos(theta2) -scale.*sin(theta2) 0; scale.*sin(theta2) scale.*cos(theta2) 0; 0 0 1]);

[scaledRotatedMoving1,RrotatedScaled1] = imwarp(K,tform1,'SmoothEdges', true);

scaledRotatedMoving1 = manage_windowing(scaledRotatedMoving1,1);

% This step is equivalent to: 
%   [scaledRotatedMoving2,RrotatedScaled2] = imwarp(moving,tform2)
% We do this to gain efficiency in computing scaledRotatedMoving2,
scaledRotatedMoving2 = rot90(scaledRotatedMoving1,2);
RrotatedScaled2 = imref2d(size(scaledRotatedMoving1),...
                          sort(-RrotatedScaled1.XWorldLimits),...
                          sort(-RrotatedScaled1.YWorldLimits));

% Form 2-D spectra associated with scaledRotatedMoving1, scaledRotatedMoving2, and fixed.
outSize = [F_tx_samples, F_tx_samples]; 
M1 = fft2(scaledRotatedMoving1,outSize(1),outSize(2));
F  = fft2(noise_frame,outSize(1),outSize(2));
M2 = fft2(scaledRotatedMoving2,outSize(1),outSize(2));

% Form the phase correlation matrix d1 for M1 correlated with F.
ABConj = F .* conj(M1);
d1 = ifft2(ABConj ./ abs(eps+ABConj),'symmetric');

% Form the phase correlation matrix d2 for M2 correlated with F.
ABConj = F .* conj(M2);
d2 = ifft2(ABConj ./ abs(eps+ABConj),'symmetric');

% Find the translation vector that aligns scaledRotatedMoving1 with fixed and
% scaledRotatedMoving2 with fixed. Choose S,theta,translation estimate that has
% the highest peak correlation in the final translation recovery step.
[vec1,peak1] = find_translation_phase_corr(d1);
[vec2,peak2] = find_translation_phase_corr(d2);

if peak1 >= peak2
    vec = vec1;
    tform = tform1;
    RrotatedScaled = RrotatedScaled1;
    peak = peak1;
else
    vec = vec2;
    tform = tform2;
    RrotatedScaled = RrotatedScaled2;
    peak = peak2;
end

% The scale/rotation operation performed prior to the final
% phase-correlation step results in a translation. The translation added
% during scaling/rotation is defined by RrotatedScaled. Form the final
% effective translation by summing the translation added during
% rotation/scale to the translation recovered in the final translation
% step.
est_shift_x  = vec(1) + (RrotatedScaled.XIntrinsicLimits(1)-RrotatedScaled.XWorldLimits(1));
est_shift_y  = vec(2) + (RrotatedScaled.YIntrinsicLimits(1)-RrotatedScaled.YWorldLimits(1));

tform.T(3,1:2) = [est_shift_x, est_shift_y];

end

%--------------------------------------------
function img = manage_windowing(img,windowing)

if windowing
    img = img .* create_blackmanWindow(size(img));
end

end

%--------------------------------------------
function h = create_blackmanWindow(windowSize)
% Define Blackman window to reduce finite image replication effects in
% frequency domain. Blackman window is recommended in (Stone, Tao,
% McGuire, Analysis of image registration noise due to rotationally
% dependent aliasing).

M = windowSize(1);
N = windowSize(2);

a0 = 7938/18608;
a1 = 9240/18608;
a2 = 1430/18608;

n = 1:N;
m = 1:M;

% Make outer product degenerate if M or N is equal to 1.
h1 = 1;
h2 = 1;
if M > 1
    h1 = a0 - a1*cos(2*pi*m / (M-1)) + a2*cos(4*pi*m / (M-1));
end
if N > 1
    h2 = a0 - a1*cos(2*pi*n / (N-1)) + a2*cos(4*pi*n / (N-1));
end

h = h1' * h2;

end
