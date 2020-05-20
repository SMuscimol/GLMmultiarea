function qf = qfregressors(dspec, toRemove, regrGroups)
%%% Returns a quadratic form that imposes smoothness on the set of
%%% regressors. More precisely, smoothness is facilitated only among
%%% regressors belonging to the same group as defined in regrGroups.
%%% dspec : design specification structure
%%% toRemove : (1xR) cell array of (1x1) cell arrays, containing strings
%%%     with the regressors to be removed.
%%% regrGroups : cell array in which each cell contains an array of indices
%%%     of regressors that are to be grouped together

if nargin < 3
    regrGroups = cell(length(dspec.covar),1);
    for r=1:length(regrGroups)
        regrGroups{r} = r;
    end
end

if nargin < 2
	toRemove = {};
end

nCovars = length(dspec.covar);
idxs = 1:nCovars;
iToR = [];
if ~isempty(toRemove)
	for r=1:length(toRemove)
		if ~isempty(toRemove{r})
			iToR = cat(1,iToR, find(strcmp({dspec.covar.label}, toRemove{r})));
		end
	end
end
idxs(iToR) = [];

qf = [ 1. ] ;
for r = 1:length(regrGroups)
    rToR = [];
    for k=1:length(iToR)
        if ~isempty(find(regrGroups{r}==iToR(k)))
            rToR = cat(1, rToR, find(regrGroups{r}==iToR(k))); 
        end
    end
    regrGroups{r}(rToR) = [];
    if ~isempty(regrGroups{r})
        if length(regrGroups{r})>1 || dspec.covar(regrGroups{r}(1)).edim>1
            qf = blkdiag(qf, qfsmooth1D(sum([dspec.covar(regrGroups{r}).edim])));
        end
        if length(regrGroups{r})==1 && dspec.covar(regrGroups{r}(1)).edim==1
            qf = blkdiag(qf, 1.);
        end
    end
end