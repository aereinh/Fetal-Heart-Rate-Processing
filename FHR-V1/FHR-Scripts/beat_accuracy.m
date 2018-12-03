%Returns a value that measures the accuracy of lineup between calculated
%beat markers and the raw signal
%
%Inputs: beats - Calculated beat markers (seconds)
%        raw_sig - Raw signal as measured by TOITU (V)
%        raw_sig_time - Time vector corresponding to raw signal with a
%           given sampling rate
%        frame_time - Range of times around each beat marker to be compared
%           to raw_sig (seconds).
%
%Outputs: acc_val - The value of accuracy of lineup between calculated beat
%           markers and raw_sig (higher=more accurate). Should be
%           interpretted as a relative value


function [acc_val] = beat_accuracy(beats,raw_sig,raw_sig_time,frame_time)

%frame_time specifies the range around each beat marker to compare to
%signal. This range should always be smaller than the difference between
%two adjacent beats
if frame_time >= min(diff(beats))
    error('frame_time should be less than minimum difference between beats');
elseif isempty(frame_time)
    frame_time = 0.1;
end

%Convert beat times and frame_time into closest number of indices that
%correspond to those times
beat_indices = round(beats./raw_sig_time(2))+1;
frame_ind_range = round(frame_time/raw_sig_time(2))+1;

%Find the mean of the absolute values of the raw_sig around the beat
%markers (within frame_range of either side of beat). Call this lineup
%value
lineup_vals = nan(1,length(beats));
for i = 1:length(beats)
    FHR_frame = abs(raw_sig(beat_indices(i)-frame_ind_range:beat_indices(i)+frame_ind_range));
    lineup_vals(i) = mean(FHR_frame);
end

%Find the mean of each beat's lineup value
acc_val = mean(lineup_vals);

end

