function [mass mass_error] = get_electron_mass
% returns electron mass in MeV/c^2
freq_shift_fit    = freq_shift();
freq_shift_slope  = freq_shift_fit.Coefficients{2, 1};

mass = 1 / freq_shift_slope;
end
