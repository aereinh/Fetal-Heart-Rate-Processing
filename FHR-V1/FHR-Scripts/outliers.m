%Simple outlier detector that bases detection on a global threshold and
%local standard deviation threshold
%
%Implementation: [percent_outliers, outlier_locs] = outliers(sig,global_thresh,nhood,loc_sd_thresh)
%
%Inputs: sig - Signal of interest
%        global_thresh - Threshold for all y values: [lower_bound,upper_bound]
%        nhood - Length of indices over which local outliers are
%           checked for; should be an odd integer > 1
%        loc_sd_thresh - Threshold for detecting local outliers; all neighborhoods 
%           with standard deviations higher than this value have outliers
%
%Outputs: percent_outliers - Percentage of total points that are labeled as
%           outliers
%         all_outlier_locs - Indexed location of all outliers

function [percent_outliers, outlier_locs] = outliers(sig,global_thresh,nhood,loc_sd_thresh)

N = length(sig);

%%%Find global outliers
if length(global_thresh) ~= 2
    global_thresh = [-Inf,Inf];
end
global_thresh = sort(global_thresh);
global_outlier_locs = sig > global_thresh(2) | sig < global_thresh(1);

%%%Find local outliers
if ~isscalar(loc_sd_thresh)
    loc_sd_thresh = Inf;
end

if isempty(nhood)
    nhood = 3;
elseif nhood <= 2 || nhood >= N
    error('nhood_length must between 3 and the length of FHR');
elseif mod(nhood,2) == 0
    nhood = nhood-1;
    disp('nhood changed to odd number: ')
    disp(nhood)
end

%Start by finding the moving standard deviation of nhood_length points
%Local neighborhoods are centered around each point with range
%nhood_length. For endpoints, available points are used

%Neighborhoods with outliers have standard deviations above cutoff
mov_sd = nan(1,N);
inc = (nhood-1)/2;

loc_outlier_locs = ones(1,N);

%Middle points
for i = inc+1:N-inc
    mov_sd(i) = std(sig(i-inc:i+inc));
    if mov_sd(i) <= loc_sd_thresh
        loc_outlier_locs(i-inc:i+inc) = 0;
    end
end

%Starting points and end points
for j = 1:inc
    mov_sd(j) = std(sig(1:j+inc));
    if mov_sd(j) <= loc_sd_thresh
        loc_outlier_locs(1:j+inc) = 0;
    end
    
    mov_sd(N-j+1) = std(sig(N-j+1-inc:end));
    if mov_sd(N-j+1) <= loc_sd_thresh
        loc_outlier_locs(N-j+1-inc:end) = 0;
    end
end

%Group local and global outliers, and find percentag
outlier_locs = find(loc_outlier_locs+global_outlier_locs > 0);
percent_outliers = length(outlier_locs)/N * 100;


end

