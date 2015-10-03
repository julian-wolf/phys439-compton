function [fit_Thomson fit_KN] = diff_cross_section (draw_figure, material)
% material can be 'BR' for brass or 'Cu' for copper
if nargin < 2
    material = ''; % default to Al, which has no label!
    if nargin < 1
        draw_figure = true;
    end
end
if nargin >= 2
    material = [material '_'];
end

[angle_offset ~] = find_zero_deg(false);

angles      = 210:5:245;
% angles      = 220:5:240;
angles_true = abs(angles - 180 - angle_offset) * pi / 180;

calib_coeffs = find_calib_coeff(false);

data_path = ['data/freq_shift_' material];

spectrum_section_210deg = ...
    normalize_spectrum(get_spectrum([data_path '210deg.chn']));
spectrum_section_215deg = ...
    normalize_spectrum(get_spectrum([data_path '215deg.chn']));
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
spectrum_section_245deg = ...
    normalize_spectrum(get_spectrum([data_path '245deg.chn']));

spectrum_section_210deg_bg = ...
    normalize_spectrum(get_spectrum([data_path '210deg_bg.chn']));
spectrum_section_215deg_bg = ...
    normalize_spectrum(get_spectrum([data_path '215deg_bg.chn']));
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
spectrum_section_245deg_bg = ...
    normalize_spectrum(get_spectrum([data_path '245deg_bg.chn']));

spectra_shift = [spectrum_section_210deg ...
                 spectrum_section_215deg ...
                 spectrum_section_220deg ...
                 spectrum_section_225deg ...
                 spectrum_section_230deg ...
                 spectrum_section_235deg ...
                 spectrum_section_240deg ...
                 spectrum_section_245deg];
             
             
spectra_section_bg = [spectrum_section_210deg_bg ...
                    spectrum_section_215deg_bg ...
                    spectrum_section_220deg_bg ...
                    spectrum_section_225deg_bg ...
                    spectrum_section_230deg_bg ...
                    spectrum_section_235deg_bg ...
                    spectrum_section_240deg_bg ...
                    spectrum_section_245deg_bg];
                
spectra_section_fg_data = double([spectra_shift] - [spectra_section_bg]);
                           
section_starts = [460 440 420 400 380 360 340 320];
% section_starts = [420 400 380 360 340];

detector_efficiency = @(E) 2.5e-12*E.^5 - 6.3e-9*E.^4 + 5.9e-6*E.^3 - ...
                           2.3e-3*E.^2 + 0.17*E + 95;

hit_counts = zeros(size(angles));

for i = 1:numel(angles)
    f = fit([section_starts(i):1024].', ...
            spectra_section_fg_data(section_starts(i):1024, i), ...
            'gauss1');
    fwhm_range = ceil(f.b1 - f.c1) : floor(f.b1 + f.c1);
    
    corrected_hits = spectra_section_fg_data(fwhm_range, i) ./ ...
        detector_efficiency( ...
            calib_coeffs(spectra_section_fg_data(fwhm_range, i)) ...
        );
    
    hit_counts(i) = sum(corrected_hits);
end

alpha0 = 1.2947; % 661.6 keV / ((511 keV / c^2) * c^2)

dOmega = 1; N = 1; I0 = 1; % TODO: use true values!

recorded_cross_sections = hit_counts / (dOmega * N * I0);

Thomson_cross_sections = 1 + (cos(angles_true)).^2;
KN_cross_sections = (1 + (cos(angles_true)).^2) ./ ...
                    ((1 + alpha0 .* (1 - cos(angles_true))).^2) .* ...
                    (   1 + ...
                        (alpha0^2 .* (1 - cos(angles_true)).^2) ./ ...
                        (   (1 + (cos(angles_true)).^2) .* ...
                            (1 + alpha0 .* (1 - cos(angles_true))) ...
                        ) ...
                    );

fit_Thomson = fitlm(Thomson_cross_sections.', recorded_cross_sections.');
fit_KN      = fitlm(KN_cross_sections.',      recorded_cross_sections.');

if draw_figure
    % Thomson figure
    figure;
    hold on
    
    sp_data = subplot('Position', [0.1 0.3 0.8 0.6]);
    plot(fit_Thomson);
    ax = gca;
    set(ax, 'TickLabelInterpreter', 'LaTeX');
    set(ax, 'xtick', []);
    % ylim([1.8 3.0]);
    ax.YTick = ax.YTick(2:end);
    ylabel('recorded cross section', 'Interpreter', 'LaTeX');
    title('');
    
    sp_resid = subplot('Position', [0.1 0.1 0.8 0.2]);
    plot(Thomson_cross_sections, fit_Thomson.Residuals.Raw, 'o');
    % plot([0.1 0.7], [0 0], '--r');
    ax = gca;
    set(ax, 'TickLabelInterpreter', 'LaTeX');
    ax.YTick = ax.YTick(2:end-1);
    ylabel('raw residuals', 'Interpreter', 'LaTeX');
    xlabel('Thomson cross section', 'Interpreter', 'LaTeX');
    title('');
    
    linkaxes([sp_data sp_resid].', 'x');
    
    % K-N figure
    figure;
    hold on
    
    sp_data = subplot('Position', [0.1 0.3 0.8 0.6]);
    plot(fit_KN);
    ax = gca;
    set(ax, 'TickLabelInterpreter', 'LaTeX');
    set(ax, 'xtick', []);
    % ylim([1.8 3.0]);
    ax.YTick = ax.YTick(2:end);
    ylabel('recorded cross section', 'Interpreter', 'LaTeX');
    title('');
    
    sp_resid = subplot('Position', [0.1 0.1 0.8 0.2]);
    plot(KN_cross_sections, fit_KN.Residuals.Raw, 'o');
    % plot([0.1 0.7], [0 0], '--r');
    ax = gca;
    set(ax, 'TickLabelInterpreter', 'LaTeX');
    ax.YTick = ax.YTick(2:end-1);
    ylabel('raw residuals', 'Interpreter', 'LaTeX');
    xlabel('Klein-Nishina cross section', 'Interpreter', 'LaTeX');
    title('');
    
    linkaxes([sp_data sp_resid].', 'x');
end
                
end
