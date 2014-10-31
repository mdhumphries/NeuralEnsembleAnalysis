NeuralEnsembleAnalysis
======================

MATLAB toolbox for analysing neural ensembles

So you've found some ensembles of neurons (or ``cell assemblies") in your population recording data - now what? This toolbox tackles this problem by laying out a set of tools for analysing neural ensembles.

It contains a collection of scripts and functions for analysing ensembles assuming:
(1) you've got a clustering of your spike-train data into ensembles (as might be obtained by: https://github.com/mdhumphries/SpikeTrainCommunitiesToolBox)
(2) you've got your spike-trains in some time-series format (i.e. either binned or "spike-density" - convolved each spike with a window function) (as might be obtained by: https://github.com/mdhumphries/SpikeTrainCommunitiesToolBox)

The focus is on classifying ensembles: this allows us to then go back to each recording and start to decompose the dynamical systems captured within
