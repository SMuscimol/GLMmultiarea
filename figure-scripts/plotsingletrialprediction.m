function [t,k] = plotsingletrialprediction(tIndex, n, expt, tIndexType, binning, mainBinsize)
%%% t index given in expt coordinates

if nargin<5
    mainBinsize = 0.1;
end
if nargin<4
    binning = 1;
end
if nargin < 3
    tIndexType = 'expt';
end

if strcmp(tIndexType,'expt')
    % change to test set coordinates %
    [t, k] = getttestfromtexpt(tIndex, n, expt);
    ti = tIndex;
elseif strcmp(tIndexType,'testSet')
    t = tIndex(1);
    k = tIndex(2);
    ti = gettexptfromttest(t,k,n,expt);
end

iS = sum(arrayfun(@(i) expt.binfun(expt.trial(ti).duration), 1:t-1)) + 1;
%iE = iS + expt.binfun(expt.trial(tIndex).duration) -1;
iE = iS + expt.binfun(expt.trial(ti).duration) -1;

yP = mainBinsize * n.yP{k};
%yP = poissrnd(yP);
yTest = full(n.yTest{k});

binnedSp = zeros(floor((iE-iS+1)/binning),1);
tRange = linspace(0.,  expt.trial(ti).duration, size(yP(iS:iE),1));
tRangeBinned = linspace(0., expt.trial(ti).duration, size(binnedSp,1));

for i=1:numel(binnedSp)
    binnedSp(i) = sum(yTest(iS + (i-1)*binning:iS + i*binning - 1));
end
binnedSp = binnedSp./(binning*mainBinsize);
%figure();
hold on
plot(tRangeBinned, binnedSp,'--', 'DisplayName','data','LineWidth',1,'Color',[0.,0.,0.]);
plot(tRange, yP(iS:iE)./mainBinsize, 'DisplayName','model','LineWidth',1,'Color',[0.3,0.3,0.3]);
%legend();
xlabel('t [s]');
ylabel('spike count [Hz]');

end