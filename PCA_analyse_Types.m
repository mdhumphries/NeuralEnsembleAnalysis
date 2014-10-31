%% decompose dynamics: 
% use ensemble types to select subpopulation of neurons in each recording,
% and do PCA of that subset
%
% Run after analysing spike-train properties and detecting ensemble types
%
% Mark Humphries 31/10/2014

clear all; close all

% load ensemble analysis data
load TestData\DataList.mat;  % list of spike files
load TestData\DataSet_ConsensusClustering_TEST_spkfcn.mat;  % convolved spike-train functions for PCA

load Analyses_Neurons_and_Groups groupdata
load Ensemble_Types AcorrTypes

nfiles = numel(DataList);
VarExplained = 0.95;


%% do PCA on selected sub-sets

for iM = 1:nfiles
    % ensemble IDs of this recording in groupdata structure
    ixGrps = find([groupdata(:).Recording] == iM);

    % ensemble-types of this subset of ensembles: 
    osctypes = zeros(numel(ixGrps),1); 
    for iG = 1:numel(ixGrps)
        osctypes(iG) = AcorrTypes(ixGrps(iG));
    end

    % select subset of groups for PCA: specify here what to use, then save
    % with appropriate names....
    PCAsubset(iM).ixGrps = ixGrps(find(osctypes == 2)); % use only groups with significant peak AND trough
    % PCAsubset(iM).ixGrps = ixGrps(find(osctypes == 3)); % bursters
    
    % now look up those neuron IDs
    PCAsubset(iM).neuronIDs = [];
    for iP = 1:numel(PCAsubset(iM).ixGrps)
        PCAsubset(iM).neuronIDs = [PCAsubset(iM).neuronIDs; groupdata(PCAsubset(iM).ixGrps(iP)).IDs];
    end
     

    % do PCA on that subset: covariance of spike-density functions
    [PCAsubset(iM).coeffs, PCAsubset(iM).scores, PCAsubset(iM).eigvalues] = princomp(spkfcn_dataset{iM}(:,PCAsubset(iM).neuronIDs));

      % compute variance explained
     PCAsubset(iM).prctVar = cumsum(PCAsubset(iM).eigvalues) / sum(PCAsubset(iM).eigvalues);
     PCAsubset(iM).nPCs = sum(PCAsubset(iM).prctVar <= VarExplained);
     
%      % reconstruct those PCs as linear sums over spike functions     
%      PCAsubset(iM).PCproj = zeros(numel(PCAsubset(iM).smples),PCAsubset(iM).nPCs);
%      for iP = 1:PCAsubset(iM).nPCs
%         smpdata = spkfcn_dataset{iM}(PCAsubset(iM).smples,PCAsubset(iM).neuronIDs); 
%         PCAsubset(iM).PCproj(:,iP) = sum(bsxfun(@times,smpdata,PCAsubset(iM).coeffs(:,iP)'),2);  % raw
%         % PCAsubset(iM).PCproj(:,iP) = sum(bsxfun(@times,bsxfun(@minus,smpdata,mean(smpdata)),PCAsubset(iM).coeffs(:,iP)'),2); % normalised to mean rate
%      end
     
     figure(iM); clf; figure(iM+100); clf
     for iP = 1:PCAsubset(iM).nPCs
         figure(iM);    title(['Recording ' num2str(iM) '; all PC projections in 95%'])
         subplot(PCAsubset(iM).nPCs,1,iP),plot(PCAsubset(iM).scores(:,iP));
         figure(iM+100);    title(['Recording ' num2str(iM) '; all factor loadings'])
         subplot(PCAsubset(iM).nPCs,1,iP),bar(PCAsubset(iM).coeffs(:,iP));      
     end
     
     if PCAsubset(iM).nPCs >= 2
         figure(iM+1000)
         plot(PCAsubset(iM).scores(:,1),PCAsubset(iM).scores(:,2)); hold on
         plot(PCAsubset(iM).scores(1,1),PCAsubset(iM).scores(1,2),'kX','MarkerSize',10)
         title(['Recording ' num2str(iM) '; % in 2 PCs = ' num2str(PCAsubset(iM).prctVar(2))])
     end

end


%% more results plots: every PC1
h = figure
hmax = figure; hold on

ixPC1 = find(arrayfun(@(x)(~isempty(x.scores)),PCAsubset,'UniformOutput', true));

ctr = 1;
for iM = ixPC1
    figure(h) 
    subplot(numel(ixPC1),1,ctr),plot(PCAsubset(iM).scores(:,1));
    title(['Recording ' num2str(iM)])
    xlabel('Sample')
    ctr = ctr+1;
    figure(hmax)
    ts = PCAsubset(iM).scores(:,1);
    plot((ts + min(ts)) ./ (max(ts) - min(ts)));  % shift to [0 1]
    xlabel('Sample')
end

%%
nNeurons = arrayfun(@(x)(numel(x.neuronIDs)),PCAsubset,'UniformOutput', true);


%%
save PCAsubsetAnalysis PCAsubset


