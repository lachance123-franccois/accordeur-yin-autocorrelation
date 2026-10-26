function [f0, conf] = pitch_autocorr(frame, fs, fmin, fmax)
frame = frame(:) - mean(frame);
N = numel(frame);
X = fft(frame, 2^nextpow2(2*N));
r = real(ifft(abs(X).^2));
r = r(1:N) / (r(1) + eps);
tmin = max(2, floor(fs / fmax));
tmax = min(N - 2, ceil(fs / fmin));
[conf, k] = max(r(tmin+1 : tmax+1));
tau = tmin + k - 1;                      
if tau > tmin && tau < tmax
    a = r(tau); b = r(tau+1); c = r(tau+2);
    den = a - 2*b + c;
    if den < 0, tau = tau + 0.5 * (a - c) / den; end
end
f0 = fs / tau;
if conf < 0.3, f0 = NaN; end
end
