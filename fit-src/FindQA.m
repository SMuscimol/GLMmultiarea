function [ActiveIdx, QuietIdx] = FindQA(HitTrials,Criteria)

%%%% Inputs:
% HitTrials is a vector including hit trial indices or hit trial numbers
% Criteria is a metrix of movement e.g. Jaw_sp of size nxp (n is the total
% number of trials and p is the number of frames which should equal to 4000
% for position or 3999 in the case of speed

%%%% Outputs:
%ActiveIdx is a binary vector of active Hit trials
%QuietIdx is a binary vector of quiet Hit trials 

%note that the size of these vectors corresponds to number of hit trials
%and not the total number of trials

DelaySpeed=nanmean(Criteria(HitTrials,500:1000),2);
BaselineSpeed=nanmean(Criteria(HitTrials,1:500),2);

[N,edges] = histcounts(BaselineSpeed,0:(max(BaselineSpeed)-min(BaselineSpeed))/100:.8*max(BaselineSpeed),'normalization','probability');
bincenters=edges(1:end-1)+(edges(2)-edges(1))/2;
[~,MAXBin]=max(N);
Thr=bincenters(MAXBin)+2*mad(BaselineSpeed);

QuietIdx=DelaySpeed<=Thr;
ActiveIdx=DelaySpeed>Thr;

end
