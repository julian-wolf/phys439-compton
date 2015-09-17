function [fit_to_data conf_int] = freq_shift (draw_final, draw_intermediate)
if nargin < 2
    draw_intermediate = false;
    if nargin < 1
        draw_final = true;
    end
end

angles      = [210 215 220 225 230 235 240 245];
angles_true = abs(angles - 180) * pi / 180;

calib_coeffs = find_calib_coeff(false);

spectrum_shift_210deg = get_spectrum('data/freq_shift_210deg.chn');
spectrum_shift_215deg = get_spectrum('data/freq_shift_215deg.chn');
spectrum_shift_220deg = get_spectrum('data/freq_shift_220deg.chn');
spectrum_shift_225deg = get_spectrum('data/freq_shift_225deg.chn');
spectrum_shift_230deg = get_spectrum('data/freq_shift_230deg.chn');
spectrum_shift_235deg = get_spectrum('data/freq_shift_235deg.chn');
spectrum_shift_240deg = get_spectrum('data/freq_shift_240deg.chn');
spectrum_shift_245deg = get_spectrum('data/freq_shift_245deg.chn');

spectrum_shift_210deg_bg = get_spectrum('data/freq_shift_210deg_bg.chn');
spectrum_shift_215deg_bg = get_spectrum('data/freq_shift_215deg_bg.chn');
spectrum_shift_220deg_bg = get_spectrum('data/freq_shift_220deg_bg.chn');
spectrum_shift_225deg_bg = get_spectrum('data/freq_shift_225deg_bg.chn');
spectrum_shift_230deg_bg = get_spectrum('data/freq_shift_230deg_bg.chn');
spectrum_shift_235deg_bg = get_spectrum('data/freq_shift_235deg_bg.chn');
spectrum_shift_240deg_bg = get_spectrum('data/freq_shift_240deg_bg.chn');
spectrum_shift_245deg_bg = get_spectrum('data/freq_shift_245deg_bg.chn');

spectra_shift = [spectrum_shift_210deg ...
                 spectrum_shift_215deg ...
                 spectrum_shift_220deg ...
                 spectrum_shift_225deg ...
                 spectrum_shift_230deg ...
                 spectrum_shift_235deg ...
                 spectrum_shift_240deg ...
                 spectrum_shift_245deg];
             
spectra_shift_bg = [spectrum_shift_210deg_bg ...
                    spectrum_shift_215deg_bg ...
                    spectrum_shift_220deg_bg ...
                    spectrum_shift_225deg_bg ...
                    spectrum_shift_230deg_bg ...
                    spectrum_shift_235deg_bg ...
                    spectrum_shift_240deg_bg ...
                    spectrum_shift_245deg_bg];
                
spectra_shift_fg_data = double([spectra_shift.data] - ...
                               [spectra_shift_bg.data]);

section_starts = [460 440 420 400 380 360 340 320];

energies    = zeros(size(angles));
frac_errors = zeros(size(angles));

for i = 1:numel(angles)
    f = fit([section_starts(i):1024].', ...
            spectra_shift_fg_data(section_starts(i):1024, i), ...
            'gauss1');
    
    if draw_intermediate
        figure;
        hold on
        
        plot(1:1024, spectra_shift_fg_data(:, i));
        plot(f);
    end
    
    
    confidences       = confint(f);
    energy_confidence = (confidences(2, 2) - confidences(1, 2)) / 2;
    
    energies(i)    = bin_to_energy(f.b1, calib_coeffs);
    frac_errors(i) = abs(bin_to_energy(energy_confidence, calib_coeffs)) / ...
                     energies(i);
end

fit_to_data = fit((1 - cos(angles_true)).', (1000 ./ (energies)).', 'poly1');
conf_int    = confint(fit_to_data);

if draw_final
    figure;
    hold on

    fill([0.001 0.599 0.599 0.001], ...
         [conf_int(1, 1) * 0.001 + conf_int(1, 2), ...
          conf_int(1, 1) * 0.599 + conf_int(1, 2), ...
          conf_int(2, 1) * 0.599 + conf_int(2, 2), ...
          conf_int(2, 1) * 0.001 + conf_int(2, 2)], ...
         'c', 'edgecolor', 'w');
    plot([0.0 0.6], fit_to_data.p1 * [0.0 0.6] + fit_to_data.p2);
    errorbar((1 - cos(angles_true)), 1000 ./ energies, ...
             frac_errors .* (1000 ./ energies), '.k');

    ylabel('inverse energy (MeV$^{-1}$)', 'Interpreter', 'LaTeX');
    xlabel('$1 - \cos \theta$',           'Interpreter', 'LaTeX');

    set(gca, 'TickLabelInterpreter', 'LaTeX');
end

end
