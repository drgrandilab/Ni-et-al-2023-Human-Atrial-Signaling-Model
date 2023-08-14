function [peaks, amplitude, latency, duration, num_pks] = get_peaks_trace(trace, min_peak_prominence,min_peak_distance)

    [peaks, latency, duration] = findpeaks(trace, 'MinPeakProminence', min_peak_prominence, 'MinPeakDistance', min_peak_distance); % two peaks closer than 100ms is treated as one peak
    latency     = latency + 1;
    
%     duration    = duration * 10;
    
    num_pks    = length(peaks);
    
    amplitude   = [];
    if num_pks >= 1
        amplitude   = peaks - min(trace);
    end
    
    num_pks    = length(peaks);

%      figure
%      plot(trace)
%      hold on
%      scatter(latency, peaks);
%      hold off
end