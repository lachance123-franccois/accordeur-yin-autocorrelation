function [name, cents, fRef, midi] = freq_to_note(f, a4)
if nargin < 2, a4 = 440; end
names = {'Do', 'Do#', 'Re', 'Re#', 'Mi', 'Fa', 'Fa#', 'Sol', 'Sol#', 'La', 'La#', 'Si'};
if isnan(f) || f <= 0
    name = '--'; cents = NaN; fRef = NaN; midi = NaN; return;
end
midiF = 69 + 12 * log2(f / a4);
midi = round(midiF);
fRef = a4 * 2^((midi - 69) / 12);
cents = 1200 * log2(f / fRef);
name = sprintf('%s%d', names{mod(midi, 12) + 1}, floor(midi / 12) - 1);
end
