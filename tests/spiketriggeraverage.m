function [minStd] = spiketriggeraverage(spikeList, signal, neuronNum, epochPreS, epochPostS, FS, potentialThreshold, plotResults)
% function [] = spiketriggeraverage(spikeList, neuronNum, signal)
% WORKS ONLY FOR THE AFFERENCES ?! Useless for afferences
%  INPUTS : - spikeList : Format: n x 3 (or 4) 
%  Each line is:    spike_time,
%                   afferent_number, (from 0 to n_neuron-1)
%
%           - signal    : can be 1D or 2D, in this case must be [nPnts*nChan]
%           - neuronNum : Neuron number
%           - epochPreS : 
%
%  OUTPUTS
%


minStd = 0;
if size(spikeList,2)==4
    if nargin<7
        error('Missing arguments: spikelist2spikerate(spikeList, nAfferents, nPnts, FS, timeBinS, potentialThreshold)');
    end
    spikeList = spikeList(spikeList(:,4)>=potentialThreshold,:);
end


nPnts   = size(signal,1);
if min(size(signal))==1
    signalSel = signal;
else
    signalSel = signal(:,neuronNum);
end
spikeListSel    = spikeList(spikeList(:,2)==neuronNum);
if isempty(spikeListSel)
    disp(['No spike for neuron ',num2str(neuronNum)]);
    return
else
    nSpikes     = length(spikeListSel);
    disp([num2str(nSpikes),' spikes for neuron ',num2str(neuronNum)]);
end
spikeTimes          = spikeListSel(:,1);
epochPreSamples     = round(epochPreS*FS);
epochPostSamples    = round(epochPostS*FS);
nPntsEpoch          = epochPreSamples+epochPostSamples+1;

triggerSig  = zeros(nSpikes,nPntsEpoch);
for i=1:nSpikes
    sig         = signalSel(max(1,spikeTimes(i)-epochPreSamples):min(nPnts,spikeTimes(i)+epochPreSamples));
    if spikeTimes(i)<=epochPreSamples
        sig     = [sig(1)*ones(epochPreSamples-spikeTimes(i)+1,1);sig];
    elseif (spikeTimes(i)+epochPostSamples)>nPnts
        sig     = [sig;sig(end)*ones(spikeTimes(i)+epochPostSamples-nPnts,1)];
    end
    triggerSig(i,:) = sig;
end
triggerSigAverage   = mean(triggerSig,1);
triggerSigNorm      = triggerSig./repmat(max(abs(triggerSig),[],2),1,nPntsEpoch);
triggerSigStd       = std(triggerSigNorm,0,1);
minStd              = min(triggerSigStd);
if minStd>1000
    disp('moko');
end
% corrBinS            = 0.03;
% corrBinSample       = round(corrBinS*FS);
% xcorrMat            = zeros(nSpikes,nPntsEpoch);
% for i=1:nSpikes
%     for t=corrBinSample:nPntsEpoch-corrBinSample
%         xcorrMat(i,t)   = max(xcorr(triggerSig(i,t:t+corrBinSample),triggerSigAverage(t:t+corrBinSample)));
%     end
% end
% xcorrMatMean = mean(xcorrMat,1);

%% Plot results
if plotResults
    figure;
    ax(1) = subplot('Position',[0.1,0.3,0.8,0.65]); hold on;
    t = linspace(-epochPreS,epochPostS,epochPostSamples+epochPreSamples+1);
    for i=1:nSpikes
        plot(t,triggerSig(i,:));
    end
    plot(t,triggerSigAverage,'color','r','linewidth',3);
    axis tight;
    plot(xlim,[0,0],'k--'); plot([0,0],ylim,'k--');
    ylabel('Amp (uV)'); set(gca,'xtick',[]);
    title(['Spike Triggered Average - neuron ',num2str(neuronNum)]);
    ax(2) = subplot('Position',[0.1,0.1,0.8,0.18]); hold on;
    plot(t,triggerSigStd);
    % plot(t,xcorrMatMean,'r');
    xlabel('Time (s)');
    linkaxes(ax,'x');
end

end
