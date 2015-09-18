function [intercept intercept_bounds] = find_zero_deg (draw_fit)
if nargin < 1
    draw_fit = false;
end

spectrum_zeroing_115deg = get_spectrum('data/Main_zeroing_115deg.chn');
spectrum_zeroing_120deg = get_spectrum('data/Main_zeroing_120deg.chn');
spectrum_zeroing_125deg = get_spectrum('data/Main_zeroing_125deg.chn');
spectrum_zeroing_130deg = get_spectrum('data/Main_zeroing_130deg.chn');
spectrum_zeroing_135deg = get_spectrum('data/Main_zeroing_135deg.chn');
spectrum_zeroing_140deg = get_spectrum('data/Main_zeroing_140deg.chn');
spectrum_zeroing_145deg = get_spectrum('data/Main_zeroing_145deg.chn');
spectrum_zeroing_150deg = get_spectrum('data/Main_zeroing_150deg.chn');

spectrum_zeroing_210deg = get_spectrum('data/Main_zeroing_210deg.chn');
spectrum_zeroing_215deg = get_spectrum('data/Main_zeroing_215deg.chn');
spectrum_zeroing_220deg = get_spectrum('data/Main_zeroing_220deg.chn');
spectrum_zeroing_225deg = get_spectrum('data/Main_zeroing_225deg.chn');
spectrum_zeroing_230deg = get_spectrum('data/Main_zeroing_230deg.chn');
spectrum_zeroing_235deg = get_spectrum('data/Main_zeroing_235deg.chn');
spectrum_zeroing_240deg = get_spectrum('data/Main_zeroing_240deg.chn');
spectrum_zeroing_245deg = get_spectrum('data/Main_zeroing_245deg.chn');

zeroing_spectra = [spectrum_zeroing_115deg spectrum_zeroing_245deg ...
                   spectrum_zeroing_120deg spectrum_zeroing_240deg ...
                   spectrum_zeroing_125deg spectrum_zeroing_235deg ...
                   spectrum_zeroing_130deg spectrum_zeroing_230deg ...
                   spectrum_zeroing_135deg spectrum_zeroing_225deg ...
                   spectrum_zeroing_140deg spectrum_zeroing_220deg ...
                   spectrum_zeroing_145deg spectrum_zeroing_215deg ...
                   spectrum_zeroing_150deg spectrum_zeroing_210deg];
                 
zeroing_spectra_normalized = zeros(numel(zeroing_spectra), 1024);
for i = 1:numel(zeroing_spectra)
    zeroing_spectra_normalized(i, :) = ...
        normalize_spectrum(zeroing_spectra(i));
end


lims65 = 360:420; % TODO
f115 = fit(lims65.', zeroing_spectra_normalized(1, lims65).', 'gauss1');
f245 = fit(lims65.', zeroing_spectra_normalized(2, lims65).', 'gauss1');
lims60 = 360:450; % TODO
f120 = fit(lims60.', zeroing_spectra_normalized(3, lims60).', 'gauss1');
f240 = fit(lims60.', zeroing_spectra_normalized(4, lims60).', 'gauss1');
lims55 = 380:480; % TODO
f125 = fit(lims55.', zeroing_spectra_normalized(5, lims55).', 'gauss1');
f235 = fit(lims55.', zeroing_spectra_normalized(6, lims55).', 'gauss1');
lims50 = 400:510;
f130 = fit(lims50.', zeroing_spectra_normalized(7, lims50).', 'gauss1');
f230 = fit(lims50.', zeroing_spectra_normalized(8, lims50).', 'gauss1');
lims45 = 420:540;
f135 = fit(lims45.', zeroing_spectra_normalized( 9, lims45).', 'gauss1');
f225 = fit(lims45.', zeroing_spectra_normalized(10, lims45).', 'gauss1');
lims40 = 440:570;
f140 = fit(lims40.', zeroing_spectra_normalized(11, lims40).', 'gauss1');
f220 = fit(lims40.', zeroing_spectra_normalized(12, lims40).', 'gauss1');
lims35 = 460:600;
f145 = fit(lims35.', zeroing_spectra_normalized(13, lims35).', 'gauss1');
f215 = fit(lims35.', zeroing_spectra_normalized(14, lims35).', 'gauss1');
lims30 = 480:630; % TODO
f150 = fit(lims30.', zeroing_spectra_normalized(13, lims30).', 'gauss1');
f210 = fit(lims30.', zeroing_spectra_normalized(14, lims30).', 'gauss1');

f_left  = fit([115 120 125 130 135 140 145].', ...
              [f115.b1 f120.b1 f125.b1 f130.b1 ...
               f135.b1 f140.b1 f145.b1].', ...
              'poly1');
f_right = fit([215 220 225 230 235 240 245].', ...
              [f215.b1 f220.b1 f225.b1 ...
               f230.b1 f235.b1 f240.b1 f245.b1].', ...
              'poly1');

conf_left  = confint(f_left);
conf_right = confint(f_right);

intercept = (f_right.p2 - f_left.p2) / (f_left.p1 - f_right.p1);

intercept_left  = (conf_right(1, 2) - conf_left(2, 2)) / ...
                  (conf_left(2, 1) - conf_right(1, 1));
intercept_right = (conf_right(2, 2) - conf_left(1, 2)) / ...
                  (conf_left(1, 1) - conf_right(2, 1));

intercept_bounds = [intercept_left intercept_right];

if draw_fit
    figure;
    hold on
    
    fill([121 239 239 121], ...
         [conf_left(1, 1) * 121 + conf_left(1, 2), ...
          conf_left(1, 1) * 239 + conf_left(1, 2), ...
          conf_left(2, 1) * 239 + conf_left(2, 2), ...
          conf_left(2, 1) * 121 + conf_left(2, 2)], ...
         'c', 'edgecolor', 'w');
    fill([121 239 239 121], ...
         [conf_right(1, 1) * 121 + conf_right(1, 2), ...
          conf_right(1, 1) * 239 + conf_right(1, 2), ...
          conf_right(2, 1) * 239 + conf_right(2, 2), ...
          conf_right(2, 1) * 121 + conf_right(2, 2)], ...
         'c', 'edgecolor', 'w');
    
    plot([115 120 125 130 135 140 145 150], ...
         [f115.b1 f120.b1 f125.b1 f130.b1 ...
          f135.b1 f140.b1 f145.b1 f150.b1], ...
         'ok');
    plot([210 215 220 225 230 235 240 245], ...
         [f210.b1 f215.b1 f220.b1 f225.b1 ...
          f230.b1 f235.b1 f240.b1 f245.b1].', ...
         'ok');
     
    plot(f_left);
    plot(f_right);
    
%     ylim([200 800]);
end

end
