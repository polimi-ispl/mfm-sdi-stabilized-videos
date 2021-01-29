addpath('..')

load('K.mat');
load('RefLogPolar_tx_K.mat');

F_tx_samples = 2^12;
F_K = fftshift(fft2(K, 2^12, 2^12));
test_F_K_logp = LogPolar_tx(F_K, [0, pi]);

if ~isequal(test_F_K_logp, F_K_logp)
    assert((sum((real(test_F_K_logp.resampledImage(:)) - ...
        real(F_K_logp.resampledImage(:))).^2)/length(test_F_K_logp.resampledImage(:)) <= eps) & ...
        (sum((imag(test_F_K_logp.resampledImage(:)) - ...
        imag(F_K_logp.resampledImage(:))).^2)/length(test_F_K_logp.resampledImage(:)) <= eps));
else
    assert(isequal(test_F_K_logp, F_K_logp))
end



