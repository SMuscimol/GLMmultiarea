function [p, stat, ps] = mywaldtest(n, dm)
%%% Runs a wald test given a neuron output structure and a dm structure.
%%% Returns the p-value and the associated statistics. Also returns the
%%% recombined p-values.

if nargin < 2 && nargout > 2
    error('must provide dm')
end

InfMatrix = inv(n.H);
variances = diag(InfMatrix);
stat = (n.wPlain)'.^2 ./ variances;
p = 1-chi2cdf(stat, 1);

if nargout == 3
    ps = combinePs(dm, p, n.wPlain(2:end));
end

end