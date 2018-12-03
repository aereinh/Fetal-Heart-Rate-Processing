%Finds the amplitude, HR, and time of the largest peak of the first
%timepoint which contains a peak within a specified range
%
%Implementation:[peak_amp,peak_HR,peak_time] = find_start_pk(spect,HR,time,start_range,quad)
%
%Inputs: spect - Spectrogram amplitude array
%        HR - Spectrogram's heart rate axis
%        time - Spectrogram's time axis
%        start_range - The range of heart rates over which starting peak is
%          searched for
%
%Outputs: peak_amps - The amplitude of the first detected peak
%         peak_HR - The heart rate that corresponding to peak_amps
%         peak_time - The time corresponding to peak_amps


function [peak_amp,peak_HR,peak_time] = find_start_pk(spect,HR,time,start_range,quad)

%Convert complex signals into real ones
if ~isreal(spect)
    spect = abs(spect);
end;

spect_size = size(spect);
lower_bound = start_range(1);
upper_bound = start_range(2);

%This while loop searches for the largest,first peak that is within the start
%range. Searches times incrementally until it finds a time that has such a
%peak.
index = [];
i = 1;
while isempty(index) && i <= spect_size(2)
    [~,locs] = findpeaks(spect(:,i),'SortStr','descend');
    test = HR(locs) >= lower_bound & HR(locs) <= upper_bound;
    index = locs(find(test,1,'first'));
    if isempty(index)
        i = i + 1;
    end
end

%If no peaks are found in the start_range, then no numbers are outputted
%If a peak is found, its position on the spectrogram is outputted
if isempty(index)
    peak_amp = nan;
    peak_HR = nan;
    peak_time = nan;
else
    peak_amp = spect(index,i);
    peak_HR = HR(index);
    peak_time = time(i);
end


%Quadratic interpolation
if ~isempty(quad) && quad~=0
    [peak_amp,peak_HR] = quad_interpolate(HR,spect(:,i),HR(index));
end
    
end