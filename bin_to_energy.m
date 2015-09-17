function [energy] = bin_to_energy(bin, fit_to_data)
fit_to_data = find_calib_coeff;
energy = fit_to_data(1) * bin + fit_to_data(2);
end
