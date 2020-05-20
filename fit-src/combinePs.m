function pout = combinePs(dm, p, w)
%%% combines the vector of p-values in a structure easy to read. It is
%%% basically an adaptation of the function 'combineWeights.m' from the
%%% neuroGLM toolbox.

dspec = dm.dspec;
binSize = dspec.expt.binSize;

if isfield(dm, 'biasCol') % undo z-score operation
    pout.bias = w(dm.biasCol);
    w(dm.biasCol) = [];
end

if isfield(dm, 'zscore') % undo z-score operation
    w = (w .* dm.zscore.sigma(:)) + dm.zscore.mu(:);
end

if isfield(dm, 'constCols') % put back the constant columns
    w2 = zeros(dm.dspec.edim, 1);
    w2(~dm.constCols) = w; % first term is bias
    w = w2;
end

if numel(w) ~= dm.dspec.edim
    error('Expecting w to be %d dimension but it''s [%d]', ...
	dspec.edim, numel(w));
end

startIdx = [1 (cumsum([dspec.covar(:).edim]) + 1)];
pout = struct();

for kCov = 1:numel(dspec.covar)
    covar = dspec.covar(kCov);
    basis = covar.basis;

    if isempty(basis)
	p_sub = p(startIdx(kCov) + (1:covar.edim) - 1);
	pout.(covar.label).tr = ((1:size(p_sub, 1))-1 + covar.offset) * binSize;
	pout.(covar.label).data = p_sub;
	continue;
    end

    assert(isstruct(basis), 'Basis structure is not a structure?');

    sdim = covar.edim / basis.edim;
    pout.(covar.label).data = zeros(size(basis.B, 1), sdim);
    for sIdx = 1:sdim
	p_sub = p(startIdx(kCov) + (1:basis.edim)-1 + basis.edim * (sIdx - 1));
	w2_sub = sum(bsxfun(@times, basis.B, p_sub(:)'), 2);
	pout.(covar.label).data(:, sIdx) = w2_sub;
    end
    pout.(covar.label).tr = ...
	(basis.tr(:, 1) + covar.offset) * binSize * ones(1, sdim);
end



end
    