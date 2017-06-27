% clear all;

%% Parameters
% eegFilepath     = 'C:\Users\deudon\Desktop\M4\_Data\EEGstim\eegstim_1chan_1patterns_100percentChan_20_200ms.mat';
% eegFilepath     = 'C:\Users\deudon\Desktop\EpiFaR\_Scripts\micMac\micmac_1_3\SAB_enc_7_tbp1_2.edf';
eegFilepath     = 'C:\Users\deudon\Desktop\M4\_Data\EEGstim\eegstim2_1patterns_50percentChan_50_50ms.edf';
fMinHz          = 20;          	% Not used if classic_v2
fMaxHz          = 200;              % Not used if classic_v2
nFilters        = 8;                % Not used if classic_v2
filterOrder     = 2;
Fe              = 2048;             % Used only if classic_v2
%- rectifier type
% rectifierType   = 'classic_v2';  	% 'classic' or 'hybrid'
rectifierType   = 'classic';        % 'classic' or 'hybrid'
% rectifierType   = 'hybrid';       % 'classic' or 'hybrid'
%- filter type
freqScale       = 'ieegEqualPower'; % 'linear' or 'erb' or 'eeg' or 'ieegEqualPower'
%- LP filter
filterOrderLP   = 2;
freqCharRatio   = 3;
% saveResults     = 1;
resultsFigDir   = 'C:\Users\deudon\Desktop\M4\Resultats\Res_main_LIFmodel';
spikeSaveDir    = 'C:\Users\deudon\Desktop\M4\Prog\matlabCochlea\mainCochlea\spikeLists';
GainiSyn        = 150;
%- STDP
N               = 64;
W               = 32;
T_i             = 8;
displayThreshold= 12;
videoOn         = 0;
%- little STDP
% N               = 4;
% W               = 4;
% T_i             = 4;
% displayThreshold= 4;


%% Load signal
if strcmpi(eegFilepath(end-2:end),'edf')
    EEG     = pop_biosig(eegFilepath);
    Fe      = EEG.srate;
    % Data sel (eeg channels)
    eegChannelInd   = cellfun(@(x)strcmp(x(1:3),'EEG'),{EEG.chanlocs.labels});
    nEegChannels    = sum(eegChannelInd);
    x               = EEG.data(eegChannelInd,:)';
    disp([num2str(nEegChannels),' EEG channels found']);
    if exist([eegFilepath(1:end-4),'_patternInd.mat'],'file')
        %- Variable is called patternInd
        load([eegFilepath(1:end-4),'_patternInd.mat']);
        x = x(1:min(length(patternInd),EEG.pnts),:);    
    end
elseif strcmpi(eegFilepath(end-2:end),'mat')
    load(eegFilepath);                  % data is named 'x'
    x               = x';
    x               = x(:,:);
    nEegChannels    = size(x,2);        % Must be only EEG channels in the data
    eegChannelInd   = 1:nEegChannels;
else
    error('Wrong input file extension');
end


%% Signal to spike conversion
tic;
[spkTimeInd, spkChanInd, vNeuron, yComp, yFilter] = signal2spikeconverter(x,Fe,nFilters,filterOrder,fMinHz,fMaxHz,freqScale,rectifierType,GainiSyn,filterOrderLP,freqCharRatio);
toc;
spikes      = vNeuron==35;


%% Save spike list
if exist([eegFilepath(1:end-4),'_patternInd.mat'],'file')
    patternId               = patternInd(spkTimeInd+1);
    patternId(patternId==0) = NaN;
    spike_list       = [double(spkTimeInd+1),double(spkChanInd+1),patternId];
else
    spike_list       = [double(spkTimeInd+1),double(spkChanInd+1),nan(size(spkTimeInd,1),1)];
end    
    
tempSep         = regexp(eegFilepath,filesep);
filename        = eegFilepath(tempSep(end)+1:end-4);
mkdir(fullfile(spikeSaveDir,filename));
spikelistPath   = fullfile(spikeSaveDir,filename,[filename,'_spike_list_0_001.mat']);
save (spikelistPath,'spike_list');


%% STDP 
tic;
output_spike_list = STDP_spike_packet_v3(spikelistPath,N,W,nFilters*nEegChannels,T_i,[],videoOn);
toc;
tic;
output_spike_list_2 = STDP_spike_packet_v4(spike_list,N,W,nFilters*nEegChannels,T_i);
toc;

%% Histogram of the different thresholds for each output spike
if ~isempty(output_spike_list)
    neuThresholds = unique(output_spike_list(:,4));
    thresholdHist = zeros(1,length(neuThresholds));
    for i=1:length(neuThresholds)
       thresholdHist(i) = sum(output_spike_list(:,4)==neuThresholds(i));
    end
    figure;
    bar(neuThresholds,thresholdHist);
else
    disp('No output spike');
end

%% STDP Results
STDP_plot_results_v3 (spike_list, output_spike_list_2, Fe, nFilters*nEegChannels, ...
    nFilters*nEegChannels, displayThreshold, size(x,1), 0.01);

return;
%% Results
% Min of std 
nNeurons        = nFilters*nEegChannels;
outputSpikeStd  = zeros(1,nNeurons);
inputSpikeStd   = zeros(1,nEegChannels);
for i=1:nNeurons
    outputSpikeStd(i) = spiketriggeraverage(output_spike_list,yFilter,i,0.2,0.2,2048,displayThreshold,0);
end
for i=1:nNeurons
    inputSpikeStd(i) = spiketriggeraverage(spike_list,yFilter,i,0.1,0.1,2048,displayThreshold,0);
end

%%
figure;
ax(1)=subplot(121);
hist(inputSpikeStd,50);
title('Input Spikes std');
ax(2)=subplot(122);
hist(outputSpikeStd,50);
title('Output Spikes std');
xlims = xlim;
axes(ax(1)); xlim(xlims);

%% Results
xDurationSec        = size(x,1)/Fe;
nPnts               = size(spikes,1);
%- Calcul of population-averaged firing rate over (10ms) time bins
% timeBinS            = 0.010;
% timeBinSample       = round(timeBinS*Fe);
% popMeanSpikeRate    = mean(spikes,2);
% nBins               = floor(nPnts/timeBinSample);
% if nBins~=nPnts/timeBinSample
%     warning('Calcul of population averaged firing rate over time bins: the end of the data will be discarded');
%     popMeanSpikeRate = popMeanSpikeRate(1:nBins*timeBinSample);
% end
% temp_popMeanSpikeRate       = reshape(popMeanSpikeRate,timeBinSample,nBins); % size TimeBinSample * nBins
% popMeanTimeMeanSpikeRate    = sum(temp_popMeanSpikeRate,1)/timeBinS;
% %- Calcul of individual firing rate 
% indiSpikeRate       = sum(spikes,1)/xDurationSec;

%- Spatio-temporal spike pattern
% figure;
% ax(1) = subplot('Position',[0.1,0.25,0.7,0.7]);
% plot(spkTimeInd/Fe, spkChanInd,'.','Color',[0.08,0.17,0.55],'MarkerSize',4);
% ylim([1,size(spikes,2)]); xlim([0,xDurationSec]);
% title('Spatio-temporal spike pattern');
% ax(2) = subplot('Position',[0.1,0.1,0.7,0.13]);
% bar(linspace(0,xDurationSec,nBins),popMeanTimeMeanSpikeRate,'stacked');
% axis tight; xlabel('Time (s)'); ylabel('Firing Rate (Hz)');
% ax(3) = subplot('Position',[0.83,0.25,0.13,0.7]);
% barh(indiSpikeRate); axis tight; xlabel('Firing Rate (Hz)');


%- Y comp
t = linspace(0,xDurationSec,length(x));
tsec2ind  	= @(x)(x-t(1))*Fe+1;
figure;
imagesc(yComp');
set(gca,'Ydir','normal');
if xDurationSec>10; tStepSec    = 10; else tStepSec=1; end;
tTickVal 	= 0:tStepSec:floor(xDurationSec);
tTickInd    = tsec2ind(tTickVal);
set(gca,'xtick',tTickInd,'xticklabel',num2str(tTickVal(:)));
ytick       = get(gca,'ytick');
% freqChar    = erbspace(fMinHz,fMaxHz,nFilters);
% yTickVal    = round(freqChar(ytick));
% set(gca,'ytick',ytick,'yticklabel',num2str(yTickVal(:)));
xlabel('time (s)'); ylabel('characteristic frequency (Hz)');
title(['Compressed signal - ',num2str(nFilters),' filters (order ',num2str(filterOrder),') - fmin: ',num2str(fMinHz),'Hz - fMax: ',num2str(fMaxHz),'Hz - FE: ',num2str(Fe),'Hz']);


%% Plot EEG Channel
% sigSel          = yComp;
sigSel          = 'filtered'; % 'filter' or 'compressed'
iChan           = 1;
if strcmpi(sigSel,'filtered')
    ySel            = yFilter(:,(iChan-1)*nFilters+1:iChan*nFilters);
elseif strcmpi(sigsel,'compressed')
    ySel            = yComp(:,(iChan-1)*nFilters+1:iChan*nFilters);
end

patternInd      = logical(patternInd);
nPatterns       = sum(diff(patternInd)==-1);
patternStarts   = find(diff(patternInd)==1);
patternEnds     = find(diff(patternInd)==-1)-1;

figure;

bx(1) = subplot(211);
plot(t,x(:,iChan)); hold on;
for i=1:nPatterns
    patternIndi = patternStarts(i):patternEnds(i);
	plot(t(patternIndi),x(patternIndi,iChan),'r');
end
axis tight
xlabel('time (s)'); ylabel('Amp (uV)');  title('Input signal');

bx(2) = subplot(212);

% Calcul the vertical space between signals
spacing         = 3*diff(ylim)/(nFilters+1);
spacingvect     = spacing.*(1:nFilters);
% Modify data to be plot in one axes
ySel        = ySel + repmat(spacingvect,length(t),1);        % Add a spacing between each signal
plot(t,ySel,'b'); hold on;
for i=1:nPatterns
    patternIndi = patternStarts(i):patternEnds(i);
	plot(t(patternIndi),ySel(patternIndi,:),'r');
end
axis tight; 
set(gca,'yticklabel',{});
xlabel('time (s)'); ylabel('Amp');  title(['Output signal - ',sigSel]);
linkaxes(bx,'x');

%% Superpose yComp for each pattern repetition
figure; 
nPatterns       = sum(diff(patternInd)==-1);
patternStarts   = find(diff(patternInd)==1);
patternEnds     = find(diff(patternInd)==-1)-1;
indPre          = round(0.02*Fe);
patternStarts   = patternStarts-indPre;
tPattern        = linspace(0,(patternEnds(1)-patternStarts(1)+1)/Fe,(patternEnds(1)-patternStarts(1)+1));

cx(1) = subplot(211); hold on;
for i=1:nPatterns
    patternIndi = patternStarts(i):patternEnds(i);
	plot(tPattern,x(patternIndi,iChan));
end
axis tight;
plot([t(indPre+1),t(indPre+1)],ylim,'k');
xlabel('time (s)'); ylabel('Amp (uV)');  title('Pattern input signal (superposed)');

cx(2)       = subplot(212); hold on;
for i=1:nPatterns
    patternIndi = patternStarts(i):patternEnds(i);
	plot(tPattern,ySel(patternIndi,:));
end
axis tight;
plot([t(indPre+1),t(indPre+1)],ylim,'k');
set(gca,'yticklabel',{});
linkaxes(cx,'x');
xlabel('time (s)'); ylabel('Amp');  title(['Pattern output signal (superposed) - ',sigSel]);