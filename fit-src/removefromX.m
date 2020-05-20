function [X, Xremoved] = removefromX(X, toRemove, dspec, removalType)
%%% Removes a certain set of regressors from X
%%% X : design matrix
%%% toRemove : (1xR) cell array of (1x1) cell arrays, containing strings
%%% with the regressors to be removed.
%%% dspec : design specification structure
%%% removalType : string. Can be 'remove' or 'shuffle'

if nargin < 4
    removalType = 'remove';
end

% find indices of regressors to be removed %
iStart = [];
iEnd = [];
for i=1:size(toRemove,2)
	if ~isempty(toRemove{i})    
    	regN = dspec.idxmap.(toRemove{i});
    	iStart(i) = sum([dspec.covar(1:(regN-1)).edim])+1;
    	iEnd(i) = iStart(i) + dspec.covar(regN).edim -1;
	end
end
iStart(find(iStart==0)) = [];
iEnd(find(iEnd==0)) = [];
idxs = [];
for i=1:length(iStart)
    idxs = horzcat(idxs, iStart(i):iEnd(i));
end

if strcmp(removalType,'remove')
    Xremoved = X(:,idxs);
    X(:,idxs) = [];
elseif strcmp(removalType,'shuffle')
    Xremoved = X(:,idxs);
    for i=idxs
        X(:,i) = X(randperm(size(X,1)),i);
    end
end

end
