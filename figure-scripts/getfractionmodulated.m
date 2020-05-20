function frac = getfractionmodulated(neurons, fitOptions, pMax)

if nargin < 3
    pMax = 0.05;
end
  
    
frac = zeros(length(fitOptions.toRemove),1);
nNeurons = length(neurons);

for r=1:length(fitOptions.toRemove)
    frac(r) = length(find(arrayfun(@(n) n.pRav(r) , neurons)<pMax))/nNeurons;
end


end