function [fit_to_data conf_int] = find_calib_coeff (draw_fit)
if nargin < 1
    draw_fit = false;
end

spectrum_133Ba_calib = get_spectrum('data/Ba133calib.chn');
spectrum_137Cs_calib = get_spectrum('data/Cs137calib.chn');
spectrum_22Na_calib  = get_spectrum('data/Na22calib.chn');

calib_spectra = [spectrum_133Ba_calib ... %
                 spectrum_133Ba_calib ... % two Ba peaks
                 spectrum_137Cs_calib ...
                 spectrum_22Na_calib];

calib_spectra_normalized = zeros(4, 1024);
for i = 1:4
    calib_spectra_normalized(i, :) = normalize_spectrum(calib_spectra(i));
end

% in KeV: Ba, Ba, Cs, Na
energies = [81 356 661.6 511];

% where to find the peaks
ranges = [70 100; 340 410; 600 800; 480 600];

bins = zeros(1, 4);
for i = 1:4
    range_crt = ranges(i, 1):ranges(i, 2);
    f = fit(range_crt.', calib_spectra_normalized(i, range_crt).', 'gauss1');
    bins(i) = f.b1;
end

fit_to_data = fit(bins.', energies.', 'poly1');
conf_int    = confint(fit_to_data);

if draw_fit
    figure;
    hold on
    
    fill([1 699 699 1], ...
         [conf_int(1, 1) *   1 + conf_int(1, 2), ...
          conf_int(1, 1) * 699 + conf_int(1, 2), ...
          conf_int(2, 1) * 699 + conf_int(2, 2), ...
          conf_int(2, 1) *   1 + conf_int(2, 2)], ...
         'c', 'edgecolor', 'w');
    plot(fit_to_data);
    plot(bins, energies, 'ok'); % error bars can't be seen either way
    
    ylabel('bin number',   'Interpreter', 'LaTeX');
    xlabel('energy (keV)', 'Interpreter', 'LaTeX');

    xlim([0 700]);
    set(gca, 'TickLabelInterpreter', 'LaTeX');
end

end
