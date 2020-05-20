function [iS, iE] = getXcovaridxs(covar, dspec)
%%% Get the start and end indices of a covariate for the design matrix X.
%%% Need the design specification structure dspec

regN = dspec.idxmap.(covar); % regressor number 
iS = sum([dspec.covar(1:(regN-1)).edim])+1;
iE = iS + dspec.covar(regN).edim -1;

end
