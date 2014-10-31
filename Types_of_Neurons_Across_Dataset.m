%%% script to analyse the clusterings of the data-set: clusters by neuron
% Mark Humphries 31/10/2014
load Analyses_Neurons_and_Groups

M = 5; % marker size
flag = '3'; % plotting flag

dists = {'sqEuclidean'};  % options for all-eigenvector consensus
rpts = 100;

minSpks = 100; % threshold: minimum number of spikes to be included in analysis

%% FILTER: only use neurons with sufficient spikes for reliable metrics in first place...

nISIs = arrayfun(@(x) numel(x.isis), neurondata);

ixKept = find(nISIs > minSpks);  % at least this many spikes
nKept = numel(ixKept);

%% clustering of data: neurons

% models are: [normal, log-normal, gamma, uniform, bivariate normal, bimodal gamma]
[nmodels tmp] = size([neurondata(ixKept).pAICs_isi]);

% (1) spike-train metrics (fit-space)
dataset_spikes = [[neurondata(ixKept).pAICs_isi]' [neurondata(ixKept).pAICs_cv2s]'];
   
[r dataD] = size(dataset_spikes);

% distance matrix
Y = pdist(dataset_spikes,'euclidean'); 
DxyData = squareform(Y);

% modularity clustering: using similarity
SxyData.Spikes = exp(-DxyData.^2);  % exponential conversion to similarity
SxyData.Spikes(eye(r)==1) = 0;      % diagonal must contain zeros...

tic
[Cneuron.Spikes,Qneuron.Spikes,Ccon_neuron.Spikes,Qc_neuron.Spikes,N_neuron,Qvec_neuron] = allevsplitConTransitive(SxyData.Spikes,dists,rpts);
toc


nGrpsCon = max(Ccon_neuron.Spikes);
for iG = 1:nGrpsCon
        % group sizes
        grpSize(iG) = sum(Ccon_neuron.Spikes == iG);
end

% [srtSize,ixSize] = sort(grpSize,'descend')
% dataset_spikes(Ccon_neuron.Spikes == ixSize(1),:) % check out in size order


%% save stuff

clear G_dataset G_data Gcon_dataset DataTable Sxy_dataset
save('Neuron_Types','Cneuron','Ccon_neuron','Qneuron','Qc_neuron','SxyData','dataset_spikes',...
    'dists','rpts','minSpks') 
