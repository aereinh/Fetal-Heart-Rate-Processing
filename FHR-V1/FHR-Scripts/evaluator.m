%Simple outlier detector that bases detection on if heart rate
%physiologically infeasible or heart rate jumps a lot. Returns location of
%outliers and frequency (%)
%
%Implementation: [percent_outliers, all_outlier_locs] = evaluator(IHR,time,HR_thresh,diff_thresh)
%
%Inputs: IHR - Calculated instantaneous heart rates at certain timepoints
%           (bpm)
%        HR_thresh - Raw threshold for HR (bpm)
%        diff_thresh - Threshold for change of HR over one second timeframe
%
%Outputs: percent_outliers - Number of outliers evaluated over total number
%           of points *100
%         all_outlier_locs - Indexed location of all evaluated outliers

function [percent_outliers, outlier_locs] = evaluator(IHR,time,HR_thresh,diff_thresh)

N = length(IHR);

%Begin by finding all values and corresponding indices of outliers
%beyond physiological range inputted by HR_thresh. Default is 100 to 200
if isempty(HR_thresh)
    HR_thresh = [100 200];
end
raw_outlier_locs = IHR > HR_thresh(end) | IHR < HR_thresh(1);

%Now look for outliers based on difference threshold. Default is Inf
if isempty(diff_thresh)
    diff_thresh = Inf;
end


%Find value in the first 5 seconds that is closest to overall median of IHR
[~,start] = min(abs(IHR(time<5)-median(IHR(time<5))));


%Assume start is not an outlier. Loop through up to start to find if a point differs
%from the last start by more than diff_thresh
outlier_locs = raw_outlier_locs;
outlier_locs(start) = 0;
if start >= 2
    for s = 1:start-1
        outlier_locs(s) = abs(IHR(s)-IHR(start)) >= (time(2)-time(1))*(sqrt(abs(s-start)))*diff_thresh;
    end
end

%Go through each point. A point is called an outlier if it differs from an
%adjacent non-outlier by more than diff_thresh. If a point is already
%classified as an outlier, base comparison on previous non-outlier
%NOTE: This method is very crude, and it could benefit a lot from
%more satisfactory techniques for non-stationary anomaly detection
for i = start+1:length(IHR)-1
    if outlier_locs(i-1) == 0 && abs(IHR(i)-IHR(i-1)) >= (time(2)-time(1))*diff_thresh
        outlier_locs(i) = 1;
    elseif outlier_locs(i-1) == 1
        ref_index = find(outlier_locs(1:i-1)==0,10,'last');
        check_diff = abs(IHR(i)-IHR(ref_index(end))) >= (time(2)-time(1))*sqrt(i-ref_index)*diff_thresh;
        if check_diff
            outlier_locs(i) = 1;
        end
    end
end

outlier_locs = find(outlier_locs==1);
percent_outliers = (length(outlier_locs)/N) * 100;

end