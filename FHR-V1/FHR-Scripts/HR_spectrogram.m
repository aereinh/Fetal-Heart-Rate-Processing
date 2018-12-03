%Uses czt to find cepstrogram (spectrogram of spectrogram of original signal)
%HR_spectrogram peaks represent fundamental frequency "candidates" (i.e. heart rates) corresponding
%to each windowed timepoint of original signal
%
%Implementation: [HR_amp,HR] = HR_spectrogram(z_sig,fz,window,HR_range)
%
%Inputs: z_sig - Array of amplitudes representing spectrogram of raw signal
%        fz - Frequency vector of raw signal's spectrogram
%        HR_range - Range of heart rates (HR = 60*fund_freq) for HR
%           spectrogram
%
%Outputs: HR_amp - Array of amplitudes corresponding to certain
%           times (columns) and heart rates (rows)
%         HR - Vector of heart rates ranging HR_range
%
%NOTE: Variables with q in front of name differentiates between spectrum of
%raw signal (no q) and the spectrum of the spectrum (with q)


function [HR_amp,HR] = HR_spectrogram(z_sig,fz,window,HR_range)

%Define z_trans to be real
if isreal(z_sig)
    z_trans = z_sig;
else
    z_trans = abs(z_sig);
end

sig_size = size(z_trans);

n = sig_size(2);
q_m = sig_size(1);

if strcmp(window,'hamming')
    window = hamming(q_m);
elseif strcmp(window,'blackman')
    window = blackman(q_m);
elseif strcmp(window,'hanning')
    window = hanning(q_m);
elseif strcmp(window,'kaiser')
    bta = input('Enter a value for beta parameter: ');
    if isempty(bta)
        bta = 0;
    end
    window = kaiser(q_m,bta);
elseif strcmp(window,'rect') || strcmp(window,'rectangle')
    window = ones(q_m,1);
elseif length(window) ~= q_m || isempty(window)
    window = ones(q_m,1);
end

%Remove DC component of z_trans to reduce excessive leakage. Window
for i = 1:n
    z_trans(:,i) = z_trans(:,i) - mean(z_trans(:,i));
    z_trans(:,i) = z_trans(:,i).*window;
end



%Define parameters of czt. Frequency range is based on HR_range paramter (reciprocal),
%qfs is based on fz (input spectrum frequency range) and czt size, qm is the size of each 
%czt bin (defined earlier)
size_fz = size(fz);

if fz(end)-fz(1) == 0
    error('Frequency range of input spectrum cannot be 0');
    
elseif size_fz(1) == 1 && size_fz(2) == 1
    q_fs = q_m/round(fz);
    
else
    q_fs = q_m/abs(fz(end)-fz(1));
end

%Default HR_range is 50 to 250 and corresponds to frequencies from 0.24 to
%1.2
if isempty(HR_range)
    q_f1 = 0.24;
    q_f2 = 1.2;
    
elseif length(HR_range)>=2
    q_f1 = 60/HR_range(end);
    q_f2 = 60/HR_range(1);
    
else
    error('Must specify a range of heart rates. If left empty, the default range is 50-250bpm');
end


%Perform optimized czt for each column vector of z_trans
q_czt = optimal_czt(z_trans,q_m,[q_f1,q_f2],q_fs);
q_fz = linspace(q_f1,q_f2,q_m);


%Invert the frequency axis of harmonic spectrum spectrum to get
%wavelength values. The values with largest amplitude correspond to most present distances
%between two sinusoidal components of signal spectrum. One of these
%corresponds to the distance between harmonics (i.e. fundamental
%frequency).
f_fund = 1./q_fz;
HR = 60*f_fund';
temp_HR_amp = abs(q_czt);

%Lastly, normalize the spectrogram to highest peak
HR_amp = nan(size(temp_HR_amp));
for k = 1:n
    HR_amp(:,k) = temp_HR_amp(:,k)/max(temp_HR_amp(:,k));
end

end