function testLL = gettestLL(yT, XT, w)
%%% Returns the log-likelihood on the test set given the design matrix XT and
%%% the model coefficients w

% set the options %
opts.family = 'poissexp';
opts.familyextra = 1;
opts.getH = true;
% format test data
xTest = full(XT);
xTest = cat(2,ones(size(xTest,1),1),xTest);
% test
opts.baseline = zeros(length(full(yT)),1,class(full(yT)));
testLL = evalGlmLikelihood(full(yT), xTest, w', opts.baseline, opts.family, 1.);

end

