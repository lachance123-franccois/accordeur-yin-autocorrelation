function [f0, t] = track_pitch(x, fs, method, frameDur, hopDur, fmin, fmax)

if nargin < 4 || isempty(frameDur), frameDur = 0.046; end
if nargin < 5 || isempty(hopDur), hopDur = 0.01; end
if nargin < 6 || isempty(fmin), fmin = 70; end
if nargin < 7 || isempty(fmax), fmax = 1000; end
x = x(:);
N = round(frameDur * fs);
hop = round(hopDur * fs);
nb = floor((numel(x) - N) / hop) + 1;
f0 = nan(1, nb);
t = ((0:nb-1) * hop + N/2) / fs;
gate = 0.01 * max(abs(x));
for b = 1:nb
    fr = x((b-1)*hop + (1:N));
    if sqrt(mean(fr.^2)) < gate, continue; end
    switch method
        case 'yin', f0(b) = pitch_yin(fr, fs, fmin, fmax);
        case 'autocorr', f0(b) = pitch_autocorr(fr, fs, fmin, fmax);
        otherwise, error('track_pitch:method', 'Methode inconnue : %s', method);
    end
end
end
