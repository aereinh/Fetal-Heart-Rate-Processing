%Removes global and local outliers and fits FHR with a cubic spline curve
%with a specified resolution
%
%Implementation: [FHR_sp,time_sp,outlier_locs] = HR_spline(FHR,time,HR_thresh,nhood,loc_sd_thresh,spline_res)
%
%Example: [FHR_sp,time_sp,outlier_locs] = HR_spline(FHR,time,[120 180],7,6,100)
%
%Inputs: FHR - Fetal heart rate extracted from HR_finder (bpm)
%        time - A vector of all timepoints used in signal (seconds)
%        nhood - An odd number of indices greater than 1 that defines the 
%           size of the local regions
%        loc_sd_thresh - Cutoff standard deviation value for finding local
%           outliers in regions defined by nhood
%        spline_res - The ratio of resolution of the spline fit to the
%           original curve
%
%Outputs: FHR_sp - Spline interpolated FHR curve, specified by outlier
%           constraints and spline resolution
%         time_sp - Time vector corresponding to FHR_sp
%         outlier_locs - Indices of outliers corresponding to given
%           constraints

function [FHR_sp,time_sp,outlier_locs] = HR_spline(FHR,time,HR_thresh,nhood,loc_sd_thresh,spline_res)

%Set default values of outlier removal and spline fitting parameters
if length(HR_thresh) ~= 2
    HR_thresh = [-Inf,Inf];
    disp('HR_thresh changed to default: HR_thresh = [-Inf,Inf]')
else
    HR_thresh = sort(HR_thresh);
end

if ~isscalar(nhood) || nhood ~= round(nhood) || nhood <= 2
    nhood = 3;
    disp('nhood changed to default: nhood = 3')
elseif nhood > length(FHR)
    error('nhood must be less than the size of FHR');
elseif mod(nhood,2) == 0
    nhood = nhood-1;
    disp('nhood changed to odd number: ')
    disp(nhood)
end

if ~isscalar(loc_sd_thresh) || loc_sd_thresh < 0
    loc_sd_thresh = Inf;
    disp('loc_sd_thresh changed to default: loc_sd_thresh = Inf')
end

if ~isscalar(spline_res) || spline_res < 1
    spline_res = 1;
    disp('spline_res changed to default: spline_res = 1')
end

%Remove outliers with outliers function if necessary based on specifications
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

end

