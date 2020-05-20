function Rsq = getRfromXy(y, X, w, yFull, mainBinsize)
% Computes the deviance-based R-sqared given a certain binned spike train
% y, the predictor matrix X, the fitted weights w. yFull is the binned
% spike train in the training set, and it is used to compute the baseline
% firing rate of the constant Poisson model. If it is not provided, the
% firing rate during the testing set is used. Mainbinsize is the size of
% each time bin used for the fit.

if nargin < 5
    mainBinsize = 0.02;
end

if nargin < 4
    yFull = y;
end

% log-likelihood %
ll = gettestLL(y, X, w);
% log-likelihood of the saturated model %
llF = evalMaxGlmLikelihood(y, 'poissexp',1);

% log-likelihood of the homogeneous model %
totNspikes = sum(full(yFull));
totTime = size(full(yFull),1)*mainBinsize ; 
homRate = totNspikes/totTime * mainBinsize;
xHom = ones(size(X,1),1);
opts.family = 'poissexp';
opts.baseline = zeros(length(full(y)),1,class(full(y)));
ll0 = evalGlmLikelihood(full(y), xHom, log(homRate), opts.baseline, opts.family, 1.);

% deviance-based R-squared %
Rsq = (ll0 - ll)/(ll0 - llF);

end