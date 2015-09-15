function [data_out] = normalize_spectrum (spectrum_in)
data_out = double(spectrum_in.data) / spectrum_in.etime;
end
