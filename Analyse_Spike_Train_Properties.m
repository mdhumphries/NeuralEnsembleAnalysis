%%% script to analyse the neuron and ensemble-level firing statistics
%
% Key outputs:
% neurondata: an array of structs, one per neuron in the whole data-set,
%             each struct containing a set of fields recording firing statistics
% groupdata: an array of structs, one per ensemble detected in the whole data-set,
%             each struct containing a set of fields recording firing statistics
%
% Key fields of both struct arrays:
%   .coeffs_isi: a cell array, one cell per fitted model, each cell a vector of coefficient values for
%                that model's maximum likelihood fit to the ISI distribution
%   .coeffs_cv2s: a cell array, one cell per fitted model, each cell a vector of coefficient values for
%                that model's maximum likelihood fit to the CV2 distribution
%   .pBICs_isis: the P(model) array for model-fits to that ISI distribution
%   .pBICs_cv2s: the P(model) array for model-fits to that CV2 distribution
%   .blnOscN: binary flag (0,1); 1 indicates a significant dip in the auto-correlogram   
%   .blnOscP: binary flag (0,1); 1 indicates a significant peak in the auto-correlogram   
%
% Mark Humphries 16/6/2014
clear all; close all

addpath Functions\  % add list of local functions to path

datapath = 'TestData/';

load([datapath 'DataSet_ConsensusClustering_TEST.mat']) % set of clustered ensembles

load([datapath 'DataList.mat']); % list of spike-train files

% select here which clustering to use: replace these fields with your data
% if using a different clustering technique
datatype = 'Con'; % 'Max'  % analyse either consensus clustering or maximum Q clusterings of the datasets
switch datatype;  
    case 'Con'   % Consensus clustering
        G_data = Gcon_dataset;  
        Ngrps = DataTable(:,7);
        Q = DataTable(:,8);
        nIDs = DataTable(:,1);
        Gbinsizes = DataTable(:,4);
    
    case 'Max'  % maximum Q clustering
        G_data = G_dataset;  
        Ngrps = DataTable(:,5);
        Q = DataTable(:,6);
        nIDs = DataTable(:,1);
        Gbinsizes = DataTable(:,4);
end

M = 5; % marker size
flag = '3'; % plotting flag

% fits to ISI and CV2 distributions
fits = [1 2 3 7 8 9]; % normal, log-normal, gamma, uniform, bivariate normal, bimodal gamma

% oscillation parameters
binsize = 1; % 1 s bins for auto-correlogram
maxlag = 20; % 20 seconds of bins...
Qt = 0.001; % quantisation step for spike-times
npermutes = 100; % number of shuffled trains for auto-correlograms

Nthresh = 2; % min. SD consecutive in trough to count as significant
Pthresh = 3; % min SD consecutive in peak to count as significant


%% analyse each neuron in clustering! 
neurondata = struct('Recording',[],'Group',[],'ID',[],'isis',[],'cv2s',[],...
                        'coeffs_isi',[],'pAICs_isi',[],'pBICs_isi',[],'nlogL_isi',[],'coeffs_cv2s',[],'pAICs_cv2s',[],'pBICs_cv2s',[],'nlogL_cv2s',[]);

groupdata = struct('Recording',[],'IDs',[],...
                        'coeffs_isi',[],'pAICs_isi',[],'pBICs_isi',[],'nlogL_isi',[],'coeffs_cv2s',[],'pAICs_cv2s',[],'pBICs_cv2s',[],'nlogL_cv2s',[]);
                    
GroupList = []; 
Sgrp = cell(numel(G_data),1); Sin = cell(numel(G_data),1); Sout = cell(numel(G_data),1);

tic
NeuronCtr = 0;
GroupCtr = 0;

for iD = 1:numel(G_data)  % cycle over all data-sets
    iD
    
    
    % make list of all groups [recording ID; group ID]: can then use this
    % to look up which groups are in which recordings after the
    % clustering.
    GroupList = [GroupList; iD*ones(Ngrps(iD),1) (1:Ngrps(iD))'];
    
    % load spike-trains
    load([datapath '\' DataList(iD).spikes]);
     
    % do sorting by similarity - assumes similarity matrix is available:
    % comment out if it isn't
    [newG,Sgrp{iD},Sin{iD},Sout{iD}] = sortbysimilarity(G_data{iD}.grps,Sxy_dataset{iD});
    
    % check for weak group membership
    [tempI,ixIn] = sort(Sin{iD},1,'descend'); [tempO,ixOut] = sort(Sout{iD}); % top-value is strongest membership
    for iN = 1:numel(tempI)
        ranks(iN,:) = [find(ixIn == iN) find(ixOut == iN)];
        totalrank(iN) = sum(ranks(iN,:));
    end
    [srtRank,ixRank] = sort(totalrank);
    
    T = [DataTable(iD,2) DataTable(iD,3)]; % time of each recording

    
%     hall = plot_clusters(spks,G_data{iD}.grps,Ngrps(iD),T,flag);  % plot in group order
%     title(['Data-set ' DataList(iD).spikes ' in group ID order (bottom-to-top)'])
    
    % first just get ISI stats for whole recording
%     IDs = unique(spks(:,1)); nIDs = numel(IDs);
%     allisis = [];
%     for iN = 1:nIDs
%         currix = find(spks(:,1) == IDs(iN));
%         ts = spks(currix,2); % spike-times of this train
%         ts = ts(ts >= T(1) & ts <= T(2)); % within the analysed period
%         allisis = [allisis; diff(ts)]; 
%     end


    % then do per-group analyses....
    for iGrps = 1:Ngrps(iD)
        GroupCtr = GroupCtr+1;
        
        thisix = find(G_data{iD}.grps(:,2) == iGrps); % array indexes of this group
        thisgrp = G_data{iD}.grps(thisix,1); % IDs of trains in this group
        nGrpIDs = numel(thisgrp);

        groupdata(GroupCtr).Recording = iD;
        groupdata(GroupCtr).IDs = thisgrp;
        
        allisis = []; allcv2 = [];         
        
        % per neuron
        for i = 1:nGrpIDs
            i
            NeuronCtr = NeuronCtr + 1;
            neurondata(NeuronCtr).Recording = iD;  % the source recording
            neurondata(NeuronCtr).Group = iGrps;   % its assigned group
            neurondata(NeuronCtr).ID = thisgrp(i); % its ID number in the recording
            
            currix = find(spks(:,1) == thisgrp(i));
            ts = spks(currix,2); % spike-times of this train
            ts = ts(ts >= T(1) & ts(ts <= T(2))); % within the analysed period
            
            neurondata(NeuronCtr).isis = diff(ts); 
            allisis = [allisis; ts(1:end-1) neurondata(NeuronCtr).isis]; % store for ensemble-level analysis
                        
            if numel(neurondata(NeuronCtr).isis) > 2
                 % MLE fits to isis                 
                [neurondata(NeuronCtr).coeffs_isi,AICs,BICs,neurondata(NeuronCtr).pAICs_isi,neurondata(NeuronCtr).pBICs_isi,neurondata(NeuronCtr).nlogL_isi] =  ...
                    fitMLEdistribution(neurondata(NeuronCtr).isis,fits);

                % CV2 of ISIs
                neurondata(NeuronCtr).cv2s = 2 * abs(neurondata(NeuronCtr).isis(1:end-1) - neurondata(NeuronCtr).isis(2:end)) ./ (neurondata(NeuronCtr).isis(1:end-1)+neurondata(NeuronCtr).isis(2:end));
                CV2_ISI(i) = mean(neurondata(NeuronCtr).cv2s);
                allcv2 = [allcv2; ts(2:end-1) neurondata(NeuronCtr).cv2s];  % store for ensemble-level analysis
                
                % MLE fits to cv2s
                cv2s_cens = neurondata(NeuronCtr).cv2s; ixCens = find(cv2s_cens <=0); cv2s_cens(ixCens) = [];
                    
                [neurondata(NeuronCtr).coeffs_cv2s,AICs,BICs,neurondata(NeuronCtr).pAICs_cv2s,neurondata(NeuronCtr).pBICs_cv2s,neurondata(NeuronCtr).nlogL_cv2s] ...
                    = fitMLEdistribution(cv2s_cens,fits);
                
                % oscillatory activity

                % auto-correlogram
                [Acorr,bins] = spike_xcorr(ts,ts,binsize,T,maxlag);
        
                % use permutation test for means for each bin
                AcorrP = [];
                for iP = 1:npermutes
                    try
                        spktsP = shuffle_intervals([ones(numel(ts),1) ts],T,Qt); 
                        tmp = spike_xcorr(spktsP(:,2),spktsP(:,2),binsize,T,maxlag);
                        AcorrP = [AcorrP; tmp];
                    catch
                        AcorrP = [];  % this ensure osc strength = 0 for this train
                        break % from this for loop
                    end
                end

                % compute means and SD
                mAcorrP = mean(AcorrP); 
                sdAcorrP = std(AcorrP); 
                ZAcorr = (Acorr-mAcorrP)./sdAcorrP; % z-transform for comparison between neurons
                ZAcorr(ceil(numel(bins)/2)) = 0;    % remove centre bin
                
                %%% consecutive bins crossing threshold
                ixsN = find(ZAcorr < -Nthresh); ixsP = find(ZAcorr > Pthresh); 
                if any(diff(ixsN)==1) 
                    neurondata(NeuronCtr).blnOscN = 1; 
                else
                    neurondata(NeuronCtr).blnOscN = 0; 
                end
                if any(diff(ixsP)==1) 
                    neurondata(NeuronCtr).blnOscP = 1; 
                else
                    neurondata(NeuronCtr).blnOscP = 0;
                end
                
            else
                % total silence....
                % probably need to fill in some blanks here...    
                
            end  % end check for enough ISIs...
        end     % end loop over spike-trains
        
        
        % fit distribution models at ensemble level: pooled ISIs and pooled
        % CV2s across all member spike-trains
        [groupdata(GroupCtr).coeffs_isi,AICs,BICs,groupdata(GroupCtr).pAICs_isi,groupdata(GroupCtr).pBICs_isi,groupdata(GroupCtr).nlogL_isi] =  ...
            fitMLEdistribution(allisis(:,2),fits);
               
        cens_cv2 = allcv2; cens_cv2(allcv2(:,2) == 0,:) = []; 
        [groupdata(GroupCtr).coeffs_cv2s,AICs,BICs,groupdata(GroupCtr).pAICs_cv2s,groupdata(GroupCtr).pBICs_cv2s,groupdata(GroupCtr).nlogL_cv2s] = ...
            fitMLEdistribution(cens_cv2(:,2),fits);      
        

        % whole ensemble auto-correlogram
        allts = sort(allisis(:,1));
        [AcorrAll,binsAll] = spike_xcorr(allts,allts,binsize,T,maxlag);
        % use permutation test for means for each bin
        AcorrP = [];
        for iP = 1:npermutes
            try
                spktsP = shuffle_intervals([ones(numel(allts),1) allts],T,Qt); 
                tmp = spike_xcorr(spktsP(:,2),spktsP(:,2),binsize,T,maxlag);
                AcorrP = [AcorrP; tmp];
            catch
                AcorrP = [];  % this ensure osc strength = 0 for this train
                break % from this for loop
            end
        end

        % compute means and SD
        mAcorrP = mean(AcorrP); 
        sdAcorrP = std(AcorrP); 
        ZAcorr = (AcorrAll-mAcorrP)./sdAcorrP; % z-transform for comparison between neurons
        ZAcorr(ceil(numel(bins)/2)) = 0;    % remove centre bin
   
        
        %%% consecutive bins crossing threshold
        ixsN = find(ZAcorr < -Nthresh); ixsP = find(ZAcorr > Pthresh); 
        groupdata(GroupCtr).thisN = any(diff(ixsN)==1);
        groupdata(GroupCtr).thisP = any(diff(ixsP)==1); 

    end % clustering loop

end % data-set loop
toc

%% SAVE ANALYSES
clear Sxy_dataset Cxy_spread_dataset Gcon_dataset G_dataset NetworkCxy_dataset stateDt_dataset binlessopts
save(['Analyses_Neurons_and_Groups'],'neurondata','groupdata','GroupList');

