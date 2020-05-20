# GLMmultiarea
Code to run GLM analysis of spiking data. Built on top of neuroGLM (https://github.com/pillowlab/neuroGLM) and also relying on the toolbox fitglmqp (https://www.mathworks.com/matlabcentral/fileexchange/31661-fit-glm-with-quadratic-penalty).

## Instructions ##
In principle, it is sufficient to execute the file fit-src/main.m, after having set the paths to the git directory, data directory, data input name and data output name (see main.m for details).
However, this runs the fit for all neurons, area, and mice group (Novice, Expert).
This takes very long (presumably days) and a lot of memory. Be sure you know what you are about to start!
It is easier to handle breaking the loop over groups and areas, and run only one combination at once.
