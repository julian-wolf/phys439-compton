function [density] = flux_density (draw_figure)
if nargin < 1
    draw_figure = false;
end

spectrum_4bricks = get_spectrum('data/flux_density_4bricks.chn');
spectrum_3bricks = get_spectrum('data/flux_density_3bricks.chn');
spectrum_2bricks = get_spectrum('data/flux_density_2bricks.chn');
spectrum_1bricks = get_spectrum('data/flux_density_1bricks.chn');
spectrum_1brick_different = ...
                   get_spectrum('data/flux_density_1brick_thinner.chn');

spectra_bricks = [spectrum_1brick_different spectrum_1bricks ...
                  spectrum_2bricks spectrum_3bricks spectrum_4bricks];
              
thicknesses = [6.7 9.3 16.0 22.7 29.4];

fit_range = 600:1024;

hit_counts = zeros(size(thicknesses));

for i = 1:5
    hits = double(spectra_bricks(i).data(fit_range));
    f = fit(fit_range.', hits, 'gauss1');
    
    fwhm_range = ceil(f.b1 - f.c1) : floor(f.b1 + f.c1);
    hit_counts(i) = sum(spectra_bricks(i).data(fwhm_range));
end

flux_fit = fit(thicknesses.', hit_counts.', 'exp1');
density  = flux_fit(0);

if draw_figure
    figure;
    hold on
    
    xlim([-2 30]);
    plot(flux_fit);
    plot(thicknesses, hit_counts, '.');
    plot(0, density, 'o');
end

end
