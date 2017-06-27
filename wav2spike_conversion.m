function [spike_list, x, Fe] = wav2spike_conversion ...
    (wavFilepath, fMinHz, fMaxHz, nFilters, filterOrder, rectifierType, bankFilterType, ...
    saveResults, resultsFigDir, GainiSyn, showFigure)
%[spike_list, x, Fe] = WAV2SPIKE_CONVERSION ...
%    (wavFilepath, fMinHz, fMaxHz, nFilters, filterOrder, rectifierType,
%    bankFilterType, saveResults, resultsFigDir, GainiSyn, showFigure)
% Converts an audio signal into a spike list, using a simple cochlea
% model. 

% To update cochlea : 
% mex -O cochlea_LIFmodel.c
% mex -O cochlea_LIFmodel_v1_5.c
% mex -O cochlea_LIFmodel_v2.c

% clear all

if nargin<11; showFigure=0; end;

%% Parameters (Script call)
if nargin==0
wavFilepath     = 'C:\Users\deudon\Desktop\M4\_Data\audioStim\speech\Tim_JAST.wav';
fMinHz          = 20;
fMaxHz          = 20000;
nFilters        = 300;
filterOrder 	= 4;    % Only multiple of 2 order
%- rectifier type
rectifierType   = 'classic';    % 'classic' or 'hybrid'
% rectifierType   = 'hybrid';     % 'classic' or 'hybrid'
%- filter type
% bankFilterType  = 'linear'; % 'linear' or 'erb'
bankFilterType  = 'erb'; % 'linear' or 'erb'
saveResults     = 1;
resultsFigDir   = 'C:\Users\deudon\Desktop\M4\_Results\Speech\linear';
GainiSyn        = 120;
end

%% Open wav file
[x,Fe]          = audioread(wavFilepath,'double');
if size(x,2)~=1;
    x = x(:,1);
end
if Fe<=2*fMaxHz; error('fMaxHz must be inferior to half the sampling frequency'); end;
xDurationSec    = length(x)/Fe;
disp(['Sampling frequency: ',num2str(Fe),' Hz.']);
disp(['Signal duration: ',num2str(xDurationSec),' sec.']);
disp(['min(x) : ',num2str(min(x)),char(10),'Max(x) : ',num2str(max(x))]);


%% Normalization of audio file
a       = 2/(max(x)-min(x));
b       = 1-a*max(x);
xNorm   = a*x+b;
x       = xNorm;
disp('Normalizing to [-1,1]');


%% Filter Design
tic;
if strcmpi(bankFilterType,'erb')
    [b,a,freqChar,ERB]  = cochlea_erbfilterbank(Fe,fMinHz,fMaxHz,nFilters,filterOrder);
elseif strcmpi(bankFilterType,'linear')
    [b,a,freqChar]      = cochlea_linearfilterbank(Fe,fMinHz,fMaxHz,nFilters,filterOrder);
else
    disp('You must specify bankFilterType (erb or linear)');
end
filterDesignTime    = toc;

%% Bank filtering with rectifier
tic;
if strcmp(rectifierType,'classic')
    [bLPfilt,aLPfilt]   = butter(1,2*65/Fe,'low');
    [Yfilt,Yrect]       = cochlea_bankfilterrectifier_classic(b,a,bLPfilt,aLPfilt,x);
elseif strcmp(rectifierType,'hybrid')
    %- Characteristic Period (for the rectifier)
    periodChar          = round(Fe./freqChar);
    periodChar          = periodChar(:);
    [Yfilt,Yrect]       = cochlea_bankfilterrectifier_hybrid(b,a,x,periodChar);
else
    error('You must specify rectifierType (erb or linear)');
end

%- Compression
Ycomp = Yrect.^(1/3);

filterAndRectifierTime = toc;

%% LIF neuron model

iSyn        = GainiSyn*Ycomp;
[vNeuron, spkTimeInd, spkChanInd] = cochlea_LIFmodel(iSyn, Fe);
spkTimeInd  = double(spkTimeInd);
spkChanInd  = double(spkChanInd);
spikes      = vNeuron==35;
% afferentNum     = spikes.*repmat(1:nFilters,size(spikes,1),1);
% tAfferentNum    = afferentNum';
% afferentNumVect = nonzeros(tAfferentNum(:));
% spikeTime       = spikes.*repmat((1:size(spikes,1))',1,nFilters);
% tSpikeTime      = spikeTime';
% spikeTimeVect   = nonzeros(tSpikeTime(:));

%- Try to open the mat file containing the pattern indices, if it exist
if exist([wavFilepath(1:end-4),'.mat'])
    %- variable is called patternInd
    load([wavFilepath(1:end-4),'.mat']);
    patternId               = nan(size(spkTimeInd,1),1);
    patternId               = patternInd(spkTimeInd+1);
    patternId(patternId==0) = NaN;
    spike_list       = [double(spkTimeInd+1),double(spkChanInd+1),patternId];
else
    spike_list       = [double(spkTimeInd+1),double(spkChanInd+1),nan(size(spkTimeInd,1),1)];
end

tempSep         = regexp(wavFilepath,filesep);
filename        = wavFilepath(tempSep(end)+1:end-4);
mkdir(fullfile('C:\Users\deudon\Desktop\M4\Prog\matlabCochlea\mainCochlea\spikeLists',filename));
save (fullfile('C:\Users\deudon\Desktop\M4\Prog\matlabCochlea\mainCochlea\spikeLists',filename,[filename,'_spike_list_0_001.mat']),'spike_list');

LIFTime = toc;
disp([num2str(length(spkTimeInd)),' spikes']);


%% Time results
disp([char(10),'Computation time ---------'])
disp(['Filter design : ',num2str(filterDesignTime),' s']);
disp(['Filtering with rectifier and compression : ',num2str(filterAndRectifierTime),' s']);
disp(['LIF model : ',num2str(LIFTime),' s']);
disp(['Total time : ',num2str(filterDesignTime+filterAndRectifierTime+LIFTime),'s']);

% return;

%% Results
if ~exist('saveResults','var'); saveResults=0; end;
if saveResults
    resultsFigDir = fullfile(resultsFigDir,filename);
    inc = 2;
    while exist(resultsFigDir,'file')
        if inc==2; 
            resultsFigDir=[resultsFigDir,'_2'];
        else
            resultsFigDir = regexprep(resultsFigDir,'_\d+$',['_',num2str(inc)]);
        end
        inc = inc+1;
    end
    mkdir(resultsFigDir);
end
    

t = linspace(0,xDurationSec,length(x));
tsec2ind  	= @(x)(x-t(1))*Fe+1;

if showFigure
    %%
    iChan = 51;
    figure;
    ax(1) = subplot(411);
    plot(t,x); title('Raw input signal');
    axis tight;
    ax(2) = subplot(412);
    plot(t,Yfilt(:,iChan)); title(['Filtered signal - channel ',num2str(iChan),' : ',num2str(round(freqChar(iChan))),' Hz']);
    ax(3) = subplot(413);
    plot(t,Yrect(:,iChan)); 
    title(['Filtered/rectified signal - channel ',num2str(iChan),' : ',num2str(round(freqChar(iChan))),' Hz']);
    ax(4) = subplot(414);
    plot(t,vNeuron(:,iChan)); hold on;
    title(['Membrane potential - channel ',num2str(iChan),' : ',num2str(round(freqChar(iChan))),' Hz']);
    xlabel('time (s)'); ylabel('Neuron potential (mV)');
    plot(xlim,[20,20],'g--');
    ylim([-5,40]);
    linkaxes(ax,'x');



    %% - Filtered and Rectified signal
    figure;
    imagesc(Ycomp');
    set(gca,'Ydir','normal');
    tStepSec    = 0.2;
    tTickVal 	= 0:tStepSec:floor(xDurationSec);
    tTickInd    = tsec2ind(tTickVal);
    set(gca,'xtick',tTickInd,'xticklabel',num2str(tTickVal(:)));
    ytick       = get(gca,'ytick');
    yTickVal    = round(freqChar(ytick));
    set(gca,'ytick',ytick,'yticklabel',num2str(yTickVal(:)));
    xlabel('time (s)'); ylabel('characteristic frequency (Hz)');
    title(['Compressed signal - ',num2str(nFilters),' filters (order ',num2str(filterOrder),') - fmin: ',num2str(fMinHz),'Hz - fMax: ',num2str(fMaxHz),'Hz - FE: ',num2str(Fe),'Hz']);
    if saveResults; 
    %     saveas(gcf,fullfile(resultsFigDir,'FilteredAndRectifiedSignal.fig')); 
        saveas(gcf,fullfile(resultsFigDir,'FilteredAndRectifiedSignal.png')); 
    end; 

    disp(['Total number of spikes: ',num2str(sum(sum(spikes==1)))]);

    %% Results - spike statistics
    % spikes = matrix (nPnts*nChan)
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
    ytick       = get(gca,'ytick');
    yTickVal    = round(freqChar(ytick));
    set(gca,'ytick',ytick,'yticklabel',num2str(yTickVal(:)));
    title('Spatio-temporal spike pattern');
    ax(2) = subplot('Position',[0.1,0.1,0.7,0.13]);
    bar(linspace(0,xDurationSec,nBins),popMeanTimeMeanSpikeRate,'stacked');
    axis tight; xlabel('Time (s)'); ylabel('Firing Rate (Hz)');
    ax(3) = subplot('Position',[0.83,0.25,0.13,0.7]);
    barh(indiSpikeRate); axis tight; xlabel('Firing Rate (Hz)');
    ytick       = get(gca,'ytick');
    yTickVal    = round(freqChar(ytick));
    set(gca,'ytick',ytick,'yticklabel',num2str(yTickVal(:)));
    linkaxes(ax(1:2),'x');
end
    
if saveResults 
    saveas(gcf,fullfile(resultsFigDir,'SpatioTemporalSpikePattern.png')); 
end

