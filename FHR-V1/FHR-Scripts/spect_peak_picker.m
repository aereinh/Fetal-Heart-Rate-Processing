%Given a heart rate spectrogram, goes through each timepoint and finds the
%heart rate and amplitude for the highest peak (with constraints)
%
%Implementation: [peak_amps,peak_HR] = spect_peak_picker(spect,HR,time,start_range,cont_constraint,quad)
%
%Inputs: spect - Spectrogram amplitude array
%        HR - Spectrogram's heart rate axis
%        time - Spectrogram's time axis
%        start_range - The range of heart rates over which the first time's peak
%          is selected
%        cont_contraint - Continuity constraint, selected peaks in a
%          certain range of the last previously detected peak
%
%Outputs: peak_amps - The amplitude of each peak heart rate for each time.
%           Times with no detected peaks are recorded as Nan
%         peak_HR - The heart rate that corresponds to the highest peak for
%           each time. Times with no detected peaks are recorded as Nan


function [peak_amps,peak_HR] = spect_peak_picker(spect,HR,time,start_range,cont_constraint,quad)

%Convert complex signals to real ones
if ~isreal(spect)
    spect_temp = spect;
    spect = abs(spect_temp);
end

%Initialize output variables to be the same size as the number of
%timepoints in the spectrogram
spect_size = size(spect);
peak_amps = nan(1,spect_size(2));
peak_HR = nan(1,spect_size(2));

%Set default start range to range of HR and default continuity constraint to no bounds
if isempty(start_range)
    start_range = [min(HR) max(HR)];
end

if isempty(cont_constraint)
    cont_constraint = Inf;
end


%Specify the first peak within the start_range
[amp1,HR1,t1] = find_start_pk(spect,HR,time,start_range,quad);

if isnan(amp1)
    error('There are no peaks in the signal over the start_range specified');
end

%Records information about starting peak in output variables
index = find(time==t1);
peak_amps(index) = amp1;
peak_HR(index) = HR1;

%Go through the following timepoints. For each one, sort the peaks and pick the
%biggest that corresponds to a frequency is within the specified cont_constraint compared to
%the previously found peak. If no peaks are in that range, that time is
%left blank (Nan) and the loop moves to the next timepoint
index = index+1;

range = start_range;
while index <= spect_size(2);
    previous_index = find(~isnan(peak_HR),1,'last');
    previous_freq = peak_HR(previous_index);    
    range(1) = previous_freq - cont_constraint;
    range(2) = previous_freq + cont_constraint;
    [ith_amp,ith_freq,ith_time] = find_start_pk(spect(:,index:end),HR,time(index:end),range,quad);
    index = find(time==ith_time);
    peak_amps(index) = ith_amp;
    peak_HR(index) = ith_freq;
    index = index + 1;
end
end