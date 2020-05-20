function nsOut = filterneurons(ns, criterion, varargin)
% ns should be a (column) struct array of neurons objects.
% criterion should be a function that takes a single neuron object and
% possibly some other inputs and returns returns a boolean that indicates
% whether the neurons satisfy the criterion

nsOut = ns(arrayfun(@(n) criterion(n,varargin{:}), ns));

end
    

