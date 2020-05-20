function duration = getTrialDuration(k, trace, videoFrameLength)
    %%% Returns the duration of trial k from the video trace 'trace'
	fullTrial = false;
	tmp = isnan(trace(k,:));
	tmp2 = diff(tmp);
	% exclude the possibility that there are some NaNs in the middle of a trial % 
	if tmp2(find(abs(tmp2)==1,1,'last'))~=-1
		firstNaN = find(abs(tmp2)==1,1,'last');
	else
		fullTrial = true;
	end

	if fullTrial
        duration = 4.; % seconds
    elseif firstNaN == 1
        duration = NaN;
    else
	duration = (firstNaN-1) * videoFrameLength ; 

    end
    
    if duration < 1 % exclude short trials
        duration = NaN;
    end
end
