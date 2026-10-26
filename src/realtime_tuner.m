function realtime_tuner(durationSec, method)
if nargin < 1 || isempty(durationSec), durationSec = 60; end
if nargin < 2 || isempty(method), method = 'yin'; end
fs = 22050;
N = 2048;                                 % ~93 ms : assez pour Mi grave (82 Hz)
rec = audiorecorder(fs, 16, 1);
record(rec);
cleanup = onCleanup(@() stop(rec));

fig = figure('Name', 'Accordeur', 'Color', 'w', 'Position', [200 200 520 360]);
ax = axes(fig, 'XLim', [-1.2 1.2], 'YLim', [-0.2 1.25]);
axis(ax, 'off'); hold(ax, 'on');
th = linspace(pi/6, 5*pi/6, 100);
plot(ax, cos(th), sin(th), 'k', 'LineWidth', 2);
for c = -50:10:50
    a = pi/2 - c/50 * pi/3;
    plot(ax, [0.92 1.0] * cos(a), [0.92 1.0] * sin(a), 'k', 'LineWidth', 1 + (c == 0));
end
needle = plot(ax, [0 0], [0 0.9], 'r', 'LineWidth', 3);
txt = text(ax, 0, -0.08, '--', 'HorizontalAlignment', 'center', 'FontSize', 28, 'FontWeight', 'bold');
sub = text(ax, 0, 1.18, '', 'HorizontalAlignment', 'center', 'FontSize', 12);

pause(0.3);
tStart = tic;
smoothC = 0;
while ishandle(fig) && toc(tStart) < durationSec
    y = getaudiodata(rec);
    if numel(y) < N, pause(0.05); continue; end
    fr = y(end-N+1:end);
    if sqrt(mean(fr.^2)) < 0.01
        set(txt, 'String', '--'); set(sub, 'String', 'jouez une note');
        pause(0.05); continue;
    end
    if strcmp(method, 'yin')
        f0 = pitch_yin(fr, fs, 70, 1000);
    else
        f0 = pitch_autocorr(fr, fs, 70, 1000);
    end
    [name, cents] = freq_to_note(f0);
    if ~isnan(cents)
        smoothC = 0.6 * smoothC + 0.4 * cents;
        a = pi/2 - max(min(smoothC, 50), -50) / 50 * pi/3;
        set(needle, 'XData', [0 0.9*cos(a)], 'YData', [0 0.9*sin(a)]);
        if abs(smoothC) < 5, col = [0 0.6 0]; else, col = [0.8 0 0]; end
        set(needle, 'Color', col);
        set(txt, 'String', name);
        set(sub, 'String', sprintf('%.1f Hz   %+.0f cents', f0, smoothC));
    end
    drawnow limitrate;
    pause(0.05);
end
end
