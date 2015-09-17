function [fit_to_data energies bins bin_errors] = find_calib_coeff
spectrum_133Ba_calib = get_spectrum('data/Ba133calib.chn');
spectrum_137Cs_calib = get_spectrum('data/Cs137calib.chn');
spectrum_22Na_calib  = get_spectrum('data/Na22calib.chn');

calib_spectra = [spectrum_133Ba_calib ... %
                 spectrum_133Ba_calib ... % two Ba peaks
                 spectrum_137Cs_calib ...
                 spectrum_22Na_calib];
             
% sources = ['133Ba'; '137Cs'; '22Na '];
             
calib_spectra_normalized = zeros(4, 1024);
for i = 1:4
    calib_spectra_normalized(i, :) = normalize_spectrum(calib_spectra(i));
end

% in KeV: Ba, Ba, Cs, Na
energies = [81 356 661.6 511];

% where to find the peaks
ranges = [70 100; 340 410; 600 800; 480 600];

bins       = zeros(1, 4);
bin_errors = zeros(1, 4);
for i = 1:4
    range = ranges(i, 1):ranges(i, 2);
    f = fit(range.', calib_spectra_normalized(i, range).', 'gauss1');
    bins(i) = f.b1;
    confidences = confint(f);
    bin_errors(i) = (confidences(2, 2) - confidences(1, 2)) / 2;
end

% figure;
% plot(bins, energies, 'o'); % error bars can't be seen either way
% hold on

fit_to_data = polyfit(bins, energies, 1);

% plot([0 700], fit_to_data(1) * [0 700] + fit_to_data(2));

end
