function [iStart, iEnd] =  getXtrialidxs(t, expt)
% Given the trial index t and the corresponding expt structure, return the
% row indices iStart and iEnd that select that trial from the corresponding
% X matrix
iStart = sum(expt.binfun([expt.trial(1:t-1).duration]))+1;
iEnd = iStart + expt.binfun(expt.trial(t).duration)-1;
end