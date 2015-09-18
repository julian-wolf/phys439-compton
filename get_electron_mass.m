function [mass mass_error] = get_electron_mass
% returns electron mass in MeV/c^2
freq_shift_fit    = freq_shift();
freq_shift_errors = confint(freq_shift_fit);
freq_shift_slope  = freq_shift_fit.p1;

slope_error            = (freq_shift_errors(2,1)-freq_shift_errors(1,1))/2;
slope_error_fractional = slope_error / freq_shift_slope;

mass       = 1 / freq_shift_slope;              % devide the given error by 2
mass_error = slope_error_fractional * mass / 2; % since it's based on 95% lims
end
