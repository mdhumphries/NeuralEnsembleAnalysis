In this directory are:
(1) Three example recordings of the Aplysia pedal ganglion during fictive locomotion.
        (i) a .mat file listing their names
        (ii) a .m file for constructing that list

Each recording's data-file contains the array "spks" that codes the spike-times of all simultaneously recorded cells in that session.

"spks" is arranged in [ID time] format (ID: neuron number; time: spike-time in seconds)

All recordings are 80 seconds long.
T=0 seconds corresponds to 1 second post-offset of the nerve stimulation that induces fictive locomotion

**************************************************
(2) A file containing an example clustering of the entire dataset into neural ensembles
Gcon_dataset is a cell array, once cell per recording. Each cell contains the ensemble assignment of each neuron in [neuronID ensembleID] format. Here the ensembles were detected using the consensus community detection algorithm from the Spike Train Community Detection Toolbox https://github.com/mdhumphries/SpikeTrainCommunitiesToolBox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(3) A file containing the convolved spike-train functions for each recording in the data-set 
Each spike convolved with a Gaussian (SD=1 s), quantised in 10 ms time-steps.





