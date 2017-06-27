
%% Parameters
eegFilepath     = 'C:\Users\deudon\Desktop\SAB\data\p27\SAB_Bipolaire\SAB 600 ms\SAB REC 1.edf';
% eegFilepath     = 'C:\Users\deudon\Desktop\SAB\data\p25\SAB_Bipolaire\SAB_ENC_4.edf';
filterOrder     = 4;
% Fe              = 2048;             % Used only if classic_v2
fMinHz          = 3;
fMaxHz          = 80;
nFilters        = 10;
%- LP filter
filterOrderLP   = 2;
freqCharRatio   = 3;
%- rectifier type
rectifierType   = 'classic_v2';     % 'classic' or 'hybrid'
freqScale       = 'ieegEqualPower'; % 'linear' or 'erb' or 'eeg' or 'ieegEqualPower'
channelSelPos   = 60:80;
% LIF
GainiSyn        = 120;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load signal
EEG     = pop_biosig(eegFilepath);
Fe      = EEG.srate;
% Data sel (eeg channels)
eegChannelInd   = cellfun(@(x)strcmp(x(1:3),'EEG'),{EEG.chanlocs.labels});
nEegChannels    = sum(eegChannelInd);
disp([num2str(nEegChannels),' EEG channels found']);

%% Channel sel
for i=1:EEG.nbchan
    disp([num2str(i),' ',EEG.chanlocs(i).labels]);
end
x   = EEG.data(channelSelPos,:)';
disp(['Keeping channel(s) : ',EEG.chanlocs(channelSelPos).labels]);

%% Normalization of the signal between -1 and 1
a 	= 2./(max(x)-min(x));
b 	= 1-a.*max(x);
x  	= x*diag(a)+repmat(b,size(x,1),1);
disp('Normalizing to [-1,1]');
nSignals = size(x,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filter Design
if strcmpi(freqScale,'erb')
    [b,a,freqChar]  = cochlea_erbfilterbank(Fe,fMinHz,fMaxHz,nFilters,filterOrder,1);
elseif strcmpi(freqScale,'linear')
    [b,a,freqChar]  = cochlea_linearfilterbank(Fe,fMinHz,fMaxHz,nFilters,filterOrder);
elseif strcmpi(freqScale,'eeg')
    [b,a,freqChar]  = cochlea_eegfilterbank(Fe,filterOrder);
    nFilters        = length(freqChar);
elseif strcmpi(freqScale,'ieegEqualPower')
    [b,a,freqChar]  = cochlea_ieegequalpowerfilterbank(Fe,fMinHz,fMaxHz,nFilters,filterOrder,1);
else
    disp('You must specify freqScale (erb or linear)');
end


%% Filtering, Rectification, Compression
Ycomp = zeros(size(x,1),nSignals*nFilters);
wb_h = waitbar(0,'Filtering, rectification and compression');
for i=1:nSignals
    waitbar(i/nSignals,wb_h);
    %% Bank filtering with rectifier
    if strcmp(rectifierType,'classic')
        [bLPfilt,aLPfilt]   = butter(1,2*65/Fe,'low');
        [~,Yrect]       = cochlea_bankfilterrectifier_classic(b,a,bLPfilt,aLPfilt,x(:,i));
    elseif strcmp(rectifierType,'classic_v2')
        aLPfilt = zeros(filterOrderLP+1,nFilters);
        bLPfilt = zeros(filterOrderLP+1,nFilters);
        for iFilt=1:nFilters
            [bLPfilt(:,iFilt), aLPfilt(:,iFilt)] = butter(filterOrderLP,(2/Fe)*freqChar(iFilt)/freqCharRatio);
        end      
        [~,Yrect]       = cochlea_bankfilterrectifier_classic_v2(b,a,bLPfilt,aLPfilt,x(:,i));
    elseif strcmp(rectifierType,'hybrid')
        %- Characteristic Period (for the rectifier)
        periodChar     	= round(Fe./freqChar);
        periodChar      = periodChar(:);
        [~,Yrect]       = cochlea_bankfilterrectifier_hybrid(b,a,x(:,i),periodChar);
    else
        error('You must specify rectifierType (erb or linear)');
    end
    %- Compression
    Yrect(Yrect<0) = 0;
    Ycomp(:,1+(i-1)*nFilters:i*nFilters) = Yrect.^(1/3);
end
close(wb_h);

%% LIF neuron model
disp('LIF model');
iSyn  = GainiSyn*Ycomp;
% [vNeuron, spkTimeInd, spkChanInd] = cochlea_LIFmodel(iSyn, Fe);
[vNeuron, spkTimeInd, spkChanInd] = cochlea_LIFmodel(iSyn, Fe);
disp([num2str(length(spkTimeInd)),' spikes']);
spikes      = vNeuron==35;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results
xDurationSec        = size(x,1)/Fe;
nPnts               = size(spikes,1);
%- Calcul of population-averaged firing rate over (10ms) time bins
timeBinS            = 0.010;
timeBinSample       = round(timeBinS*Fe);
popMeanSpikeRate    = mean(spikes,2);
nBins               = floor(nPnts/timeBinSample);
if nBins~=nPnts/timeBinSample
    warning('Calcul of population averaged firing rate over time bins: the end of the data will be discarded');
    popMeanSpikeRate = popMeanSpikeRate(1:nBins*timeBinSample);
end
temp_popMeanSpikeRate       = reshape(popMeanSpikeRate,timeBinSample,nBins); % size TimeBinSample * nBins
popMeanTimeMeanSpikeRate    = sum(temp_popMeanSpikeRate,1)/timeBinS;
%- Calcul of individual firing rate 
indiSpikeRate       = sum(spikes,1)/xDurationSec;

%- Spatio-temporal spike pattern
figure;
ax(1) = subplot('Position',[0.1,0.25,0.7,0.7]);
plot(spkTimeInd/Fe, spkChanInd,'.','Color',[0.08,0.17,0.55],'MarkerSize',4);
ylim([1,size(spikes,2)]); xlim([0,xDurationSec]);
title('Spatio-temporal spike pattern');
ax(2) = subplot('Position',[0.1,0.1,0.7,0.13]);
bar(linspace(0,xDurationSec,nBins),popMeanTimeMeanSpikeRate,'stacked');
axis tight; xlabel('Time (s)'); ylabel('Firing Rate (Hz)');
ax(3) = subplot('Position',[0.83,0.25,0.13,0.7]);
barh(indiSpikeRate); axis tight; xlabel('Firing Rate (Hz)');

%% Results - Y comp
t = linspace(0,xDurationSec,length(x));
tsec2ind  	= @(x)(x-t(1))*Fe+1;
figure;
imagesc(Ycomp');
set(gca,'Ydir','normal');
tStepSec    = 10;
tTickVal 	= 0:tStepSec:floor(xDurationSec);
tTickInd    = tsec2ind(tTickVal);
set(gca,'xtick',tTickInd,'xticklabel',num2str(tTickVal(:)));
ytick       = get(gca,'ytick');
% yTickVal    = round(freqChar(ytick));
% set(gca,'ytick',ytick,'yticklabel',num2str(yTickVal(:)));
xlabel('time (s)'); ylabel('characteristic frequency (Hz)');
title(['Compressed signal - ',num2str(nFilters),' filters (order ',num2str(filterOrder),') - FE: ',num2str(Fe),'Hz']);


