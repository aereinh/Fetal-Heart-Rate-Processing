%Converts calculated FHR into beat sinusoid
%
%Implementation: [sine_function,beats,time_shift,phase_shift,acc_val] = 
%   bpm2sines(t,IHR,raw_sig,raw_sig_time,frame_time)
%
%Inputs: t - Time vector for calculated FHR (seconds)
%        IHR - The calculated FHR for each timepoint in t (bpm)
%        raw_sig - Raw signal as measured by TOITU (V)
%        raw_sig_time - Time vector corresponding to raw signal with a
%           given sampling rate
%        frame_time - Range of times around each beat marker to be compared
%           to raw_sig (seconds)
%
%Outputs: sine_function - A complex sinusoid that represents location of
%           each heartbeat over time as described by calculated IHR
%         beats - Calculated event markers of where beats occur (seconds)
%         time_shift - The amount beats and sine_function is shifted to get
%           the most accurate lineup with raw_sig (seconds)
%         acc_val - The value of accuracy of lineup between calculated beat
%           markers and raw_sig (higher=more accurate). Should be
%           interpretted as a relative value

function [sine_function,beats,time_shift,phase_shift,acc_val] = bpm2sines(t,IHR,raw_sig,raw_sig_time,frame_time)
%The IHR at each timepoint corresponds to frequencies
freqs = IHR./60;
omega = 2*pi.*freqs;

%Integrate omega from 0 to t. Then create a complex sinusoid whose
%frequency changes over time according to IHR
phis = cumsum(omega*(t(2)-t(1)));
sine_function = exp(1j.*phis);

%Set the times corresponding to peaks of the sinusoid as the beats. These
%may need to be shifted over by some amount to maximize accuracy
[~,beat_locs] = findpeaks(real(sine_function));
beats = t(beat_locs);

%If last three inputs are specified, we can determine heuristically how
%much the beats should be shifted to maximize accuracy. Shift the sinusoid correspondingly. 
%See beat_accuracy function for more details
if ~isempty(raw_sig) && ~isempty(raw_sig_time) && ~isempty(frame_time)
    m_diff = max(diff(beats));
    if beats(1)-m_diff/2-frame_time >= 0
        count = linspace(-max(diff(beats))/2,max(diff(beats))/2,50);
    else
        count = linspace(frame_time,max(diff(beats))+frame_time,50);
    end
    
    accuracy = nan(1,length(count));
    for i = 1:length(count)
        accuracy(i) = beat_accuracy(beats+count(i),raw_sig,raw_sig_time,frame_time);
    end
    
    [acc_val,max_index] = max(accuracy);
    [~,time_shift] = quad_interpolate(count,accuracy,count(max_index));
    
    beats = beats+time_shift;
    index_shift = round(time_shift/t(2))+1;
    phase_shift = -sum(omega(1:beat_locs(1)-index_shift)*t(2)-t(1));
    sine_function = sine_function*exp(-1j*phase_shift);
    
    
%If values are not specified, don't shift beats
else
    time_shift = [];
    acc_val = [];
end

end