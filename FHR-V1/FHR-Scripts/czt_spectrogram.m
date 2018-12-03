%Chirp z transform spectrogram over certain frequency range.
%Based on MATLAB's spectrogram function.
%
%Implementation: [z_trans,fz,time] = czt_spectrogram(signal,window,noverlap,nczt,f,fs)
%
%Inputs: signal - input signal
%        window - specify type and length
%        noverlap - amount of overlap of each window (used to specify time
%          increments)
%        nczt - size of czt for each windowed region of signal
%        f - [f1,f2], frequency range of czt
%        fs - sampling rate
%
%Outputs: z_trans - matrix of amplitudes corresponding to certain times (each column)
%           and frequency (each row)
%         fz - vector of frequencies, used as one axis for spectrogram
%         time - vector of times, used as one axis for the spectrogram


function [z_trans,fz,time] = czt_spectrogram(signal,window,noverlap,nczt,f,fs)

%Window input can be a positive integer, in which case a rectangular window of
%that size is used. If the window is not positive or not an integer, it
%is not used
if isscalar(window) && (window < 1 || window ~= round(window))
    error('The inputted window must be a positive integer or a window vector');

elseif isscalar(window)
    s = window;
    window = ones(s,1);
end

%The window size must be at most the size of the signal and at least 1
if length(window) > length(signal)
    error('The window is too large for the signal');
    
elseif length(window) < 1
    error('The window must be at least size 1');
end

%noverlap must an integer between zero and length(window)-1
if noverlap ~= round(noverlap)
    error('The value of noverlap must be an integer');

elseif noverlap < 0
    error('The value of noverlap must be non-negative');
    
elseif noverlap >= length(window)
    error('The value of noverlap must be an integer less than the size of the window');
end

%Define size variables, then truncate the signal to match noverlap -- i.e. so that the total number of
%"bins" fits exactly within the signal
L_temp = length(signal);
W = length(window);
tot_bins = floor((L_temp-W)/(W-noverlap))+1;
signal = signal(1:tot_bins*(W-noverlap)+noverlap);
L = length(signal);

%Make sure signal is a row vectors
sig_size = size(signal);
if sig_size(1)~=1
    signal = signal';
end

%Split signal into an array containing all overlapping bins of the signal.
%Multiply each bin by window function.
signal_bins = zeros(W,tot_bins);
for idx = 1:tot_bins
    i_start = 1+(idx-1)*(W-noverlap);
    i_end = i_start+W-1;
    signal_bins(:,idx) = signal(i_start:i_end)'.*window;
end

    
%Define parameters of chirp z-transform

%Frequency range: f1, f2
f1 = f(1);
f2 = f(2);

%Size: m
m = nczt;
if isempty(nczt)
    m = W;

%czt runs faster if m is a power of 2, especially power of 2 less than half
%the window size
elseif any([strcmp(nczt,'optimal') strcmp(nczt,'opt') strcmp(nczt,'fast')]) == 1
    m = 2^(floor(log(W)/log(2))-1);
end

%Perform optimized czt
z_trans = optimal_czt(signal_bins,m,f,fs);
for k = 1:tot_bins
    z_trans(:,k) = z_trans(:,k)/max(z_trans(:,k));
end

fz = linspace(f1,f2,m)';
time = (linspace(0,L-W,tot_bins)/fs)';
end