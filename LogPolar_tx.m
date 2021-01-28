function I_log_polar = LogPolar_tx(I,thetaRange)
% Function inspired to 'images.internal.LogPolar'
% I_log_polar = LogPolar_tx(I,thetaRange) computes the log polar
% transformation of the image I. 
% INPUTS: 
% I = image I (gray scale)
% thetaRange = angular values for the log-polar transformation
% OUTPUT:
% I_log_polar = Log polar transformation of the image I. Contains
% information about the number of samples along rho and theta dimensions,
% the minimum and maximum values of rho and theta. 
% I_log_polar.resampledImage is the actual log-polar tx of the image I

inSize = size(I);

% Set the maximum radius to be the largest hypotenuse. When
% computing the polar FFT, we don't want to throw away
% information at the high frequency regions of the FFT.
center = 0.5 + inSize/2;
cx = center(2);
cy = center(1);

rhoMin = 1;
rhoMax = hypot(cx-0.5,cy-0.5);
% Allow numerics to degenerate reasonably if grid size is 1
% pixel along a particular dimension.
rhoMax = max(1+eps,rhoMax);

I_log_polar.numSamplesRho = round(rhoMax);

I_log_polar.thetaMin = thetaRange(1);
I_log_polar.thetaMax = thetaRange(2);

% The following equation describes the constraining
% relationship necesssary to guarantee that a pixel's nearest
% neighbors in orthogonal directions are equally spaced:
%
% rhoMin = rhoMax * exp(-2*pi*(numSamplesRho-1) / numSamplesTheta );
%
% Which can be rearranged to calculate the number of angular
% samples necessary to achieve equally spaced orthogonal
% neighbors given rhoMin,rhoMax, and numSamplesRho:
I_log_polar.numSamplesTheta = -2*pi * (I_log_polar.numSamplesRho-1) / log(rhoMin/rhoMax);

% Allow numerics to degenerate reasonably if grid size is 1
% pixel along a particular dimension.
I_log_polar.numSamplesTheta = max(1,I_log_polar.numSamplesTheta);

deltaTheta = I_log_polar.thetaMax / I_log_polar.numSamplesTheta;
theta = linspace(0,I_log_polar.thetaMax-deltaTheta,I_log_polar.numSamplesTheta);

% We sample rho such that the maximum logRho is log(rhoMax) and
% the minimum value of logRho is 1. We sample linearly in
% log-rho. It can be proven that this type of scaling is
% equivalent to an exponential of base rhoMax:
%
%    k = 1:numSamplesRho;
%    rhoBaseN = rhoMax .^ ( k ./ numSamplesRho);
%
% That is, our spacing in rho is independent of the choice of log
% base.

I_log_polar.logRhoMin = log(rhoMin);
I_log_polar.logRhoMax = max(1e-6,log(rhoMax));
logRho = linspace(I_log_polar.logRhoMin,I_log_polar.logRhoMax,I_log_polar.numSamplesRho);
rho   = exp(logRho);

% Build mesh of rho, theta
[theta,rho] = meshgrid(theta,rho);

% Find cooresponding intrinsic coordinates in input image
[X,Y] = pol2cart(theta,rho);

% Adjust X and Y to account for center of polar transformation
X = X+cx;
Y = Y+cy;

% Resample and maintain the resampled image grid in single
% precision floating point for speed unless supplied data was
% in double.
if ~isa(I,'double')
    I = single(I);
end

I_log_polar.resampledImage = images.internal.interp2d(I,X,Y,'bilinear',0);

end