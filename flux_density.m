function [density] = flux_density (draw_figure)
if nargin < 1
    draw_figure = false;
end

calib_coeffs = find_calib_coeff(false);

spectrum_4bricks = ...
    normalize_spectrum(get_spectrum('data/flux_density_4bricks.chn'));
spectrum_3bricks = ...
    normalize_spectrum(get_spectrum('data/flux_density_3bricks.chn'));
spectrum_2bricks = ...
    normalize_spectrum(get_spectrum('data/flux_density_2bricks.chn'));
spectrum_1bricks = ...
    normalize_spectrum(get_spectrum('data/flux_density_1bricks.chn'));
spectrum_1brick_different = ...
    normalize_spectrum(get_spectrum('data/flux_density_1brick_thinner.chn'));

spectra_bricks = [spectrum_1brick_different spectrum_1bricks ...
                  spectrum_2bricks spectrum_3bricks spectrum_4bricks];
              
thicknesses = [6.7 9.3 16.0 22.7 29.4] ./ 100;

fit_range = 600:1024;

detector_efficiency = @(E) (2.5e-12*E.^5 - 6.3e-9*E.^4 + 5.9e-6*E.^3 - ...
                            2.3e-3*E.^2 + 0.17*E + 95) ./ 100;

hit_counts = zeros(size(thicknesses));

for i = 1:5
    hits = double(spectra_bricks(fit_range, i));
    f = fit(fit_range.', hits, 'gauss1');
    fwhm_range = ceil(f.b1 - f.c1) : floor(f.b1 + f.c1);
    
    relevant_data  = double(spectra_bricks(fwhm_range, i));
    corrected_data = relevant_data ./ ...
                     detector_efficiency(calib_coeffs(relevant_data));
    hit_counts(i)  = sum(corrected_data);
end

flux_fit = fit(thicknesses.', hit_counts.', 'exp1');

detector_area = 0.0254^2; % as measured in the lab

density  = flux_fit(0) / detector_area;

if draw_figure
    figure;
    hold on
    
    xlim([-2 30] ./ 100);
    plot(flux_fit);
    plot(thicknesses, hit_counts, '.');
    plot(0, density, 'o');
end

end
