%Goes through HR_sig, cuts it up into regions (for faster processing), finds
%heart rate for each of those, and stiches into one FHR vector
%
%Implementation: [FHR,time,FHR_sp,time_sp,outlier_locs] = 
%HR_finder(raw_signal,region_size,spect_window,delta_time,freq_range,fs,cepst_window,HR_range,spline)
%
%Example: [FHR,time,FHR_sp,time_sp,outlier_locs] =
%HR_finder(heart_sig1,50,hamming(3.5*fs),0.5,[5,25],fs,rectwin(3.5*fs),[100,200],[120,180,8,100]);
%
%Inputs: raw_signal - Raw input signal
%        region_size - Length of regions which signal is split into to
%           perform czt_spectrogram (seconds). Too long may lead to long
%           runtime
%        Window_time - Length of windows (seconds)
%        delta_time - Spacing of time, equal to (window-noverlap)/fs (seconds)
%        fs - Sampling frequency of HR_signal (Hz)
%        spline - Gives user option to specify outlier bounds and
%           resolution for a spline interpolant FHR curve. Specify (in
%           order) global threshold, neighborhood, local sd threshold, and spline res with an
%           array, or enter 1 for prompt
%
%Outputs: FHR - The heart rate that corresponds to the highest peak of
%           HR_spectrogram for each time (bpm)
%         time - A vector of all timepoints used in signal (seconds)
%         FHR_sp - Spline interpolated FHR curve, specified by outlier
%           constraints and spline resolution
%         time_sp - Time vector corresponding to FHR_sp
%         outlier_locs - Indices of outliers corresponding to given
%           constraints



function [FHR,time,FHR_sp,time_sp,outlier_locs] = HR_finder(raw_signal,region_size,spect_window,delta_time,freq_range,fs,cepst_window,HR_range,spline)

window_time = length(spect_window)/fs;

if delta_time > window_time
    error('delta_time should be less than the window time');
end

%Calculate how many timepoints will be used. Truncate the signal to match
%exactly those number of timepoints
n_timepoints = floor(length(raw_signal)/(delta_time*fs));
raw_signal = raw_signal(1:n_timepoints*delta_time*fs);
FHR = nan(1,n_timepoints);
time = nan(1,n_timepoints);

%Take the Hilbert transform to get smoother envelope of signal
raw_signal = abs(hilbert(raw_signal));

%Define some default parameters
if isempty(region_size) || strcmp(region_size,'all')
    region_size = length(raw_signal)/fs;
end

if isempty(cepst_window)
    cepst_window = 'rect';
elseif length(cepst_window)~=length(spect_window)
    error('Windows for spectrum and cepstrum must be same size');
end

if isempty(freq_range)
    freq_range = [5,25];
elseif length(freq_range)~=2
    error('freq_range be a vector with two numbers, where the first is the smallest');
end

if isempty(HR_range)
    HR_range = [100 200];
elseif length(HR_range)~=2
    error('HR_range be a vector with two numbers, where the first is the smallest');
end

%Go through each chunk in the signal. Calculate the czt_spectrogram,
%HR_spectrogram, and highest HR's corresponding to each window of each
%timepoint in each chunk. Stitch all chunks together
i_start = 0;
i_end = i_start+region_size+window_time;
time_scale = 1;
while i_end*fs <= length(raw_signal)
    [z,fz,t] = czt_spectrogram(raw_signal(i_start*fs+1:i_end*fs),spect_window,(window_time-delta_time)*fs,'fast',freq_range,fs);
    [HR_amp,HR] = HR_spectrogram(z,fz,cepst_window,HR_range);
    [~,bpm] = spect_peak_picker(HR_amp,HR,t,[],[],1);
    FHR(time_scale:time_scale+length(t)-1) = bpm;
    if time_scale == 1
        time(1:length(t)) = t;
    else
        time(time_scale:time_scale+length(t)-1) = t+time(find(~isnan(time),1,'last'))+delta_time;
    end
    time_scale = time_scale+length(t);
    i_start = i_start+(region_size+delta_time);
    i_end = i_start+region_size+window_time; 
end

%Finish FHR off with remainder of signal
if i_start+1<length(raw_signal) && length(raw_signal(i_start*fs+1:end))>=length(spect_window)
    [z,fz,t] = czt_spectrogram(raw_signal(i_start*fs+1:end),spect_window,(window_time-delta_time)*fs,'fast',freq_range,fs);
    [HR_amp,HR] = HR_spectrogram(z,fz,cepst_window,HR_range);
    [~,bpm] = spect_peak_picker(HR_amp,HR,t,[],[],1);
    FHR(time_scale:time_scale+length(t)-1) = bpm;
    if time_scale == 1
        time(1:length(t)) = t;
    else
        time(time_scale:time_scale+length(t)-1) = t+time(find(~isnan(time),1,'last'))+delta_time;
    end
end
 
%Remove lagging nan's from FHR and time vector
nan_index = find(isnan(FHR),1,'first');
FHR(nan_index:end) = [];
time(nan_index:end) = [];

%Set default values of outlier removal and spline fitting parameters
HR_thresh = HR_range;
nhood = 3;
loc_sd_thresh = Inf;
spline_res = 1;


%If spline input is an array, the first four inputs are interpretted as:
%[lower bound, upper bound, jump constraint , spline resolution] =
%[HR_thresh,diff_thresh,spline_res]
if length(spline) >= 2
    HR_thresh = [spline(1),spline(2)];
end
if length(spline) >= 3
    nhood = spline(3);
end
if length(spline) >= 4
    loc_sd_thresh = spline(4);
end
if length(spline) >= 5
    spline_res = spline(5);
end

%Prompt user for specifications if 1 is entered
if spline == 1
    HR_thresh_temp = input('Enter the HR bounds in a 2-tuple: ');
    if length(HR_thresh_temp) == 2
        HR_thresh = HR_thresh_temp;
    end
    nhood_temp = input('Enter the neighborhood size (odd number of indices >= 3: ');
    if isscalar(nhood_temp)
        nhood = nhood_temp;
    end
    loc_sd_thresh_temp = input('Enter the local std constraint: ');
    if isscalar(loc_sd_thresh_temp)
        loc_sd_thresh = loc_sd_thresh_temp;
    end
    spline_res_temp = input('Enter the spline resolution (length of spline divided by length of original): ');
    if isscalar(spline_res_temp) && spline_res_temp >= 1
        spline_res = spline_res_temp;
    end
end

%Remove outliers with outliers function if necessary based on specifications
if spline~= 0 && (sum(HR_thresh==HR_range) > 0 || loc_sd_thresh < range(HR_range))
    [~,outlier_locs] = outliers(FHR,HR_thresh,nhood,loc_sd_thresh);
    FHR_abr = FHR;
    FHR_abr(outlier_locs) = nan;
    
    %Avoid exterpolation by removing tailing nans at beginning and end
    first_num = find(~isnan(FHR_abr),1,'first');
    if first_num > 1
        FHR_abr(1:first_num) = [];
        time(1:first_num) = [];
    end
    
    last_num = find(~isnan(FHR_abr),1,'last');
    if last_num < length(FHR_abr)
        FHR_abr(last_num:end) = [];
        time(last_num:end) = [];
    end
    
    time_sp = time(1):(time(2)-time(1))/spline_res:time(end);
    FHR_csape = csape(time,FHR_abr);
    FHR_sp = ppval(FHR_csape,time_sp);
    
else
    outlier_locs = [];
    time_sp = [];
    FHR_sp = [];
end

end