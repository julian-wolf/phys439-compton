function [energy] = bin_to_energy (bin, fit_to_data)
if nargin == 1
    fit_to_data = find_calib_coeff(false);
end

energy = fit_to_data.p1 * bin + fit_to_data.p2;

end
