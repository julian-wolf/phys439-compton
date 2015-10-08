function [norm_factor] = diff_cross_section (material)
% material can be 'BR' for brass or 'Cu' for copper
if nargin < 1
    material = ''; % default to Al, which has no label!
else
    material = [material '_'];
end

[angle_offset ~] = find_zero_deg(false);

angles      = 210:5:245;
if ~strcmp(material, '')
    angles  = 220:5:240;
end

angles_true = abs(angles - 180 - angle_offset) * pi / 180;

calib_coeffs = find_calib_coeff(false);

data_path = ['data/freq_shift_' material];

if strcmp(material, '')
    spectrum_section_210deg = ...
        normalize_spectrum(get_spectrum([data_path '210deg.chn']));
    spectrum_section_215deg = ...
        normalize_spectrum(get_spectrum([data_path '215deg.chn']));
end
spectrum_section_220deg = ...
    normalize_spectrum(get_spectrum([data_path '220deg.chn']));
spectrum_section_225deg = ...
    normalize_spectrum(get_spectrum([data_path '225deg.chn']));
spectrum_section_230deg = ...
    normalize_spectrum(get_spectrum([data_path '230deg.chn']));
spectrum_section_235deg = ...
    normalize_spectrum(get_spectrum([data_path '235deg.chn']));
spectrum_section_240deg = ...
    normalize_spectrum(get_spectrum([data_path '240deg.chn']));
if strcmp(material, '')
    spectrum_section_245deg = ...
        normalize_spectrum(get_spectrum([data_path '245deg.chn']));
end

if strcmp(material, '')
    spectrum_section_210deg_bg = ...
        normalize_spectrum(get_spectrum([data_path '210deg_bg.chn']));
    spectrum_section_215deg_bg = ...
        normalize_spectrum(get_spectrum([data_path '215deg_bg.chn']));
end
spectrum_section_220deg_bg = ...
    normalize_spectrum(get_spectrum([data_path '220deg_bg.chn']));
spectrum_section_225deg_bg = ...
    normalize_spectrum(get_spectrum([data_path '225deg_bg.chn']));
spectrum_section_230deg_bg = ...
    normalize_spectrum(get_spectrum([data_path '230deg_bg.chn']));
spectrum_section_235deg_bg = ...
    normalize_spectrum(get_spectrum([data_path '235deg_bg.chn']));
spectrum_section_240deg_bg = ...
    normalize_spectrum(get_spectrum([data_path '240deg_bg.chn']));
if strcmp(material, '')
    spectrum_section_245deg_bg = ...
        normalize_spectrum(get_spectrum([data_path '245deg_bg.chn']));
end

spectra_section = [spectrum_section_220deg ...
                   spectrum_section_225deg ...
                   spectrum_section_230deg ...
                   spectrum_section_235deg ...
                   spectrum_section_240deg];     

spectra_section_bg = [spectrum_section_220deg_bg ...
                      spectrum_section_225deg_bg ...
                      spectrum_section_230deg_bg ...
                      spectrum_section_235deg_bg ...
                      spectrum_section_240deg_bg];
                
if strcmp(material, '')
    spectra_section = [spectrum_section_210deg ...
                       spectrum_section_215deg ...
                       spectra_section ...
                       spectrum_section_245deg];
                   
    spectra_section_bg = [spectrum_section_210deg_bg ...
                          spectrum_section_215deg_bg ...
                          spectra_section_bg ...
                          spectrum_section_245deg_bg];
end

spectra_section_fg = double([spectra_section] - [spectra_section_bg]);
                           
section_starts = [460 440 420 400 380 360 340 320];
if ~strcmp(material, '')
    section_starts = [420 400 380 360 340];
end

detector_efficiency = @(E) (2.5e-12*E.^5 - 6.3e-9*E.^4 + 5.9e-6*E.^3 - ...
                            2.3e-3*E.^2 + 0.17*E + 95) ./ 100;

n_iter = numel(angles);
if ~strcmp(material, '')
    n_iter = n_iter - 1;
end

hit_counts = zeros([1 n_iter]);
count_errs = zeros([1 n_iter]);

for i = 1:n_iter
    f = fit([section_starts(i):1024].', ...
            spectra_section_fg(section_starts(i):1024, i), ...
            'gauss1');
    fwhm_range = ceil(f.b1 - f.c1) : floor(f.b1 + f.c1);
    
    left_range  = ceil(f.b1 - 4.5*f.c1) : floor(f.b1 - 2.5*f.c1);
    right_range = ceil(f.b1 + 2.5*f.c1) : floor(f.b1 + 3.0*f.c1);
    
    fit_range_bg = [left_range right_range];
    f_bg = fitlm(fit_range_bg, spectra_section_fg(fit_range_bg, i));
    
    correct_efficiency = @(counts) counts ./ ...
        detector_efficiency(calib_coeffs(counts));
    
    corrected_hits = correct_efficiency(spectra_section_fg(fwhm_range, i));
    
    [background, bg_ci] = predict(f_bg, fwhm_range.');
    
    corrected_bg = correct_efficiency(background);
    
    corrected_bg_min = correct_efficiency(bg_ci(:, 1));
    corrected_bg_max = correct_efficiency(bg_ci(:, 2));
    
    hit_counts(i) = sum(corrected_hits) - sum(corrected_bg);
    count_errs(i) = (sum(corrected_bg_max) - sum(corrected_bg_min));
    
    if true
        figure;
        hold on
        
        corrected_fit = correct_efficiency(f(1:1024));
        plot(corrected_fit);
        plot(fwhm_range, corrected_bg,     '-k' );
        plot(fwhm_range, corrected_bg_min, '--r');
        plot(fwhm_range, corrected_bg_max, '--r');
        
        figure;
        hold on
        
        plot(1:1024, spectra_section_fg(:, i));
        plot(f);
        plot(f_bg);
    end
end

alpha0 = 1.2947; % 0.19; % 661.6 keV / ((511 keV / c^2) * c^2)
dOmega = 0.165;
N      = 1.01e25;
if strcmp(material, 'BR_')
    N  = 3.18e25;
end
if strcmp(material, 'Cu_')
    N  = 2.82e25;
end
I0     = flux_density(false);
r0     = 0.2; % just a guess, really

Thomson_cross_section = @(angles) 1 + (cos(angles)).^2;
KN_cross_section      = @(angles) ...
                         (1 + (cos(angles)).^2) ./ ...
                         ((1 + alpha0 .* (1 - cos(angles))).^2) .* ...
                         (  1 + ...
                            (alpha0^2 .* (1 - cos(angles)).^2) ./ ...
                            (   (1 + (cos(angles)).^2) .* ...
                                (1 + alpha0 .* (1 - cos(angles))) ...
                            ) ...
                         );

recorded_cross_sections = hit_counts / (dOmega * N * I0) / (r0^2 / 2);
recorded_errors         = count_errs / (dOmega * N * I0) / (r0^2 / 2);

norm_fit = fitlm(recorded_cross_sections, ...
                 KN_cross_section(angles_true(1:n_iter)), 'y ~ x1 - 1')
norm_factor = norm_fit.Coefficients{1, 1};

recorded_cross_sections = recorded_cross_sections * norm_factor;
recorded_errors         = recorded_errors         * norm_factor;

all_angles = 0.0:0.01:1.9;

figure;
hold on

% sp_data = subplot('Position', [0.1 0.3 0.8 0.6]);
plot(all_angles, KN_cross_section(all_angles), '-');
% hold on
plot(all_angles, Thomson_cross_section(all_angles), '--');
errorbar(angles_true(1:n_iter), recorded_cross_sections, recorded_errors, '.');
ax = gca;
set(ax, 'TickLabelInterpreter', 'LaTeX');
% set(ax, 'xtick', []);
% ax.YTick = ax.YTick(2:end);
ylabel('recorded cross section', 'Interpreter', 'LaTeX');
title('');

% sp_resid = subplot('Position', [0.1 0.1 0.8 0.2]);
% plot(angles_true(1:n_iter), ...
%      KN_cross_section(angles_true(1:n_iter)) - ...
%      recorded_cross_sections(1:n_iter), ...
%      'o');
% hold on
% plot([0.0 1.5], [0 0], '--r');
% ax = gca;
% set(ax, 'TickLabelInterpreter', 'LaTeX');
% ax.YTick = ax.YTick(1:end-1);
% ylabel('raw residuals', 'Interpreter', 'LaTeX');
xlabel('$\theta$ (rad)', 'Interpreter', 'LaTeX');
% title('');
% 
% linkaxes([sp_data sp_resid].', 'x');
                
end
