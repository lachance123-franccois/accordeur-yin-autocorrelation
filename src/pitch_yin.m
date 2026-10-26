function [f0, aperiodicite, dprime] = pitch_yin(frame, fs, fmin, fmax, thr)
if nargin < 5 || isempty(thr), thr = 0.1; end
frame = frame(:);
N = numel(frame);
tmax = min(floor(N/2), ceil(fs / fmin));
tmin = max(2, floor(fs / fmax));
W = N - tmax;                           
d = zeros(tmax + 1, 1);
x0 = frame(1:W);
for tau = 1:tmax
    e = x0 - frame(1+tau : W+tau);
    d(tau + 1) = e' * e;
end
dprime = ones(tmax + 1, 1);
cs = cumsum(d(2:end));
dprime(2:end) = d(2:end) .* (1:tmax)' ./ (cs + eps);

idx = find(dprime(tmin+1 : end) < thr, 1);
if isempty(idx)
    m = min(dprime(tmin+1 : end));
    idx = find(dprime(tmin+1 : end) < m + 0.1, 1);
end
tau = tmin + idx - 1;
while tau + 1 <= tmax && dprime(tau + 2) < dprime(tau + 1)
    tau = tau + 1;                       
end
aperiodicite = dprime(tau + 1);
t = tau;
if tau > 1 && tau < tmax
    a = dprime(tau); b = dprime(tau + 1); c = dprime(tau + 2);
    den = a - 2*b + c;
    if den > 0, t = tau + 0.5 * (a - c) / den; end
end
f0 = fs / t;
if aperiodicite > 0.5, f0 = NaN; end
end
