%%% script to detect the types of ensembles in the data-set
%
% Key outputs:
% dataset_spikes: the matrix defining the fit-space - one row per ensemble,
%                 first set of columns the ISI P(model) vector, 
%                 the second set of columns the CV2 P(model) vector
% SxyData.Spikes: the similarity matrix for all pairs of ensembles in the
%                 fit-space
% Ccon.Spikes:    the ensemble-types identified by the consensus community
%                 detection algorithm (format: [EnsembleType#]
% Cmax.Spikes:    the ensemble-types identified by the maximum modularity (Q) clustering from the community
%                 detection algorithm (format: [EnsembleType#]
% AcorrTypes:     the ensemble-types idenitified by significant peaks or
%                 troughs in their autocorrelogram
% 
% Mark Humphries 16/6/2014
clear all; close all

% path to per-neuron and per-ensemble spike-train statistics
load Analyses_Neurons_and_Groups

M = 5; % marker size
flag = '3'; % plotting flag

dists = {'sqEuclidean'};  % options for all-eigenvector consensus
rpts = 100;

%% classify by auto-correlogram

AcorrTypes = zeros(length(groupdata),1);

% four possible outcomes: 
% (1) no significant peaks or troughs = "non-oscillatory"
AcorrTypes([groupdata.thisP] == 0 & [groupdata.thisN] == 0) = 1;  

% (2) significant peaks and troughs = "oscillator"
AcorrTypes([groupdata.thisP] == 1 & [groupdata.thisN] == 1) = 2;  

% (3) significant peaks but no significant troughs = "burster"
AcorrTypes([groupdata.thisP] == 1 & [groupdata.thisN] == 0) = 3;  

% (4) no significant peaks but significant troughs = "pauser"
AcorrTypes([groupdata.thisP] == 0 & [groupdata.thisN] == 1) = 4;  


%% spike-train metrics: fits to distributions

% concatenate vectors of model-fits: make fit-space
dataset_spikes = [[groupdata.pAICs_isi]' [groupdata.pAICs_cv2s]'];

Y = pdist(dataset_spikes,'euclidean'); % find distance between each pair of ensembles in fit-space
DxyData = squareform(Y);
nnodes = size(DxyData,1)

% modularity clustering: using similarity
SxyData.Spikes = exp(-DxyData.^2);  % exponential conversion to similarity
SxyData.Spikes(eye(nnodes)==1) = 0; % zeros on diagonal

% cluster using consensus community detection

tic
[Cmax.Spikes,Qmax.Spikes,Ccon.Spikes,Qc.Spikes,N,Qvec] = allevsplitConTransitive(SxyData.Spikes,dists,rpts);
toc

% Ccon.Spikes IDs ensemble-type of each ensemble across recordings

%% save stuff
clear G_dataset G_data Gcon_dataset DataTable Sxy_dataset ensembleRate neurondata 

save('Ensemble_Types','Cmax','Ccon','Qmax','Qc','SxyData','dataset_spikes','dists','rpts','AcorrTypes')


