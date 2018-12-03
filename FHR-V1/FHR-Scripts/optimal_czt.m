%Optimized version of czt, based on Grant D Martin's algorithm
%'Chirp Z-Transform Spectral Zoom Optimization with MATLAB®'
%http://prod.sandia.gov/techlib/access-control.cgi/2005/057084.pdf
%
%Implementation: [z_trans] = optimal_czt(x,m,f,fs)
%
%Inputs: x - Input signal
%        m - Length of czt. Type 'fast' for an efficient choice of m
%        f - Frequency range over which czt is performed (Hz)
%        fs - Sampling frequency of signal (Hz)
%
%Ouputs: z_trans - czt transform of x

function [z_trans] = optimal_czt(x,m,f,fs)

%Begin Spectral Zoom Optimization CZT
f1 = f(1);
f2 = f(2);

[k, n] = size(x);
oldk = k;
if k == 1
    x = x(:);
    [k, n] = size(x);
end

if isempty(m)
    m = length(x);
end

%Length for power of 2 fft
nfft = 2^nextpow2(k+m-1);

%Pre-multiply data
kk = ((-k+1):max(m-1,k-1)).';
kk2 = (kk.^2)./2;
ww = exp(-1i*2*pi*(f2-f1)/(m*fs)*kk2);
nn = (0:(k-1))';
aa = exp(-1i*2*pi*f1/fs.*nn);
aa = aa.*ww(k+nn);
y = x.*aa(:,ones(1,n));

%Fast convolution via FFT
fy = fft(y,nfft);
fv = fft(1./ww(1:(m-1+k)),nfft);
fy = fy.*fv(:,ones(1,n));
z_trans = ifft(fy);

%Final multiply
z_trans = z_trans(k:(k+m-1),:).*ww(k:(k+m-1),ones(1,n));

if oldk == 1
    z_trans = z_trans';
end
end