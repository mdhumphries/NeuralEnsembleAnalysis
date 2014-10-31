function ctrlspkts = shuffle_intervals(spkts,T,dt)

% SHUFFLE_INTERVALS shuffle inter-event times of neural time-series
% C = SHUFFLE_INTERVALS(S,T,DT) returns the two-column time-series vector C [ID,time], containing shuffled
% events for each neuron. Vector S contains the original
% time-series; T is the two-element vector giving [start end] times of the
% recording in seconds; DT is the quantisation step of the original spike-times (in seconds).
%
% Possible that the function could get stuck in an infinite loop in
% pathological cases; so returns with error message after 20 attempts
% 
% Mark Humphries 6/10/2012

ctrlspkts = [];
cellIDs = unique(spkts(:,1));

for i = 1:numel(cellIDs)    
    Train = spkts(spkts(:,1)==cellIDs(i),2);        % onset times in seconds
    IEI = diff(Train);                              % inter-event intervals for train 
    IStrt = Train(1)-T(1); IEnd = T(2) - Train(end);  % time before first event; time after last
    ntime_free = IStrt + IEnd;                      % total amount of free time before first and after last event
    nintervals = numel(IEI);
    blnOK = 0; nattempts = 0;
    while ~blnOK
        new_ix = randperm(nintervals);                    % randomly shuffle indices of intervals in train
        rndStrt = T(1) + ceil(rand * ntime_free / dt) * dt;      % randomly chosen start time, quantised to original spike-time resolution
        shufTs = [rndStrt; rndStrt+cumsum(IEI(new_ix))];  % starting from randomly chosen start time, the times of new events in shuffled train 1,
        if max(shufTs) > T(2) warning('last event in shuffled train occurs too late');   % paranoia - should rarely happen, except with pathological cases (too few spikes)
            nattempts = nattempts + 1;
            if nattempts > 20
                error('Cannot satisfactorily shuffle this spike-train');
            end
        else blnOK = 1;
        end
    end
    ctrlspkts = [ctrlspkts; zeros(nintervals+1,1)+cellIDs(i) shufTs];   % add event indices, add to shuffled onset times vector
end
