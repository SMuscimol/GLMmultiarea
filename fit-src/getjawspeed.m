function jawspeed = getjawspeed(trace , binSize, trial, expt, fitParams)
%%% Return the jaw speed variable, computed by averaging the mean of the
%%% absolute value on a bin of size binSize

% remove NaNs
trace(find(isnan(trace))) = 0.;

jawspeed = zeros(expt.binfun(trial.duration),1);
regBinSize = round(binSize./fitParams.mainBinsize);

for i=1:floor(expt.binfun( trial.duration)/regBinSize) 
    jawspeed((i-1)*regBinSize+1:i*regBinSize) = mean(abs(...
        trace((i-1)*round(binSize/fitParams.videoFrameLength)+1:min(size(trace,1),i*round(binSize/fitParams.videoFrameLength)),1)));    
end

if isnan(jawspeed(end,1))
    jawspeed(end,1)=0.;
end

end
