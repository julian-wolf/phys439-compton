function [spectrum] = get_spectrum (filename)
[fd, msg] = fopen(filename,'r');
if fd < 0
    msg = [msg, ' error returned by attempt to open "', filename, '"'];
    error('get_spectrum:openError',msg);
end

spectrum = get_spectrum_chn(fd);
fclose(fd);

end

function [spectrum] = get_spectrum_chn (fd)
type = fread(fd,1,'int16=>int16');
if type ~= -1
    error('get_spectrum_chn:headerError','Error reading header');
end

unit_number    = fread(fd, 1, 'uint16=>unit16');
segment_number = fread(fd, 1, 'uint16=>unit16');
ascii_seconds  = fread(fd, 2, 'uint8=>char');
real_time_20ms = fread(fd, 1, 'uint32=>unit32');
live_time_20ms = fread(fd, 1, 'uint32=>unit32');
start_date     = fread(fd, 8, 'uint8=>char');
start_time     = fread(fd, 4, 'uint8=>char');
channel_offset = fread(fd, 1, 'uint16=>unit16');

spectrum.numbchanspm   = fread(fd, 1, 'uint16=>uint16');
spectrum.data          = fread(fd, spectrum.numbchanspm, 'uint32=>uint32');
spectrum.pca_date      = [start_date(3:5, 1)', ' ', start_date(1:2, 1)', ...
                          ' ', sprintf('%02d', 19+start_date(8, 1)-'0'), ...
                          start_date(6:7, 1)'];
spectrum.acqstart_date = spectrum.pca_date;

hours   =  sscanf(start_time(1:2, 1), '%d');
minutes = sscanf(start_time(3:4, 1), '%d');
seconds = sscanf(ascii_seconds(1:2, 1), '%d');

spectrum.acqstart_time = sprintf('%02d:%02d:%02d %cm',...
                               1+mod(hours-1, 12),  minutes, ...
                               seconds, 'a' + (hours > 12) .* ('p'-'a'));
spectrum.acquire_start = 24 * 3600 * datenum(spectrum.acqstart_date) + ...
                         3600 * hours + 60 * minutes + seconds;
spectrum.acquire_stop  = spectrum.acquire_start + 0.02 * real_time_20ms;
spectrum.etime         = 0.02 * live_time_20ms;
end
