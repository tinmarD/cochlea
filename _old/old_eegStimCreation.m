%%% EEG Stim creation
warning('The script will create a signal WITH discontinuities');

rng('shuffle');

%% Parameters
inputSignalPath     = 'C:\Users\deudon\Desktop\EpiFaR\_Scripts\micMac\micmac_1_3\SAB_enc_7_tbp1_2.edf';
patternStart        = 40;       % sec 
patternDuration     = 0.2;      % sec
noiseDurationInt    = [0.2 0.5];            % sec
patternChanProp     = 100;      % percentage
nRepets             = 20;
patternColors       = {'c','g','r','b','y'};
nPatterns           = length(patternDuration);


%% Load signal
EEG = pop_biosig(inputSignalPath);
%- Remove non-eeg channels
eegChannelInd   = cellfun(@(x)strcmp(x(1:3),'EEG'),{EEG.chanlocs.labels});
eegChannelPos   = nonzeros(eegChannelInd.*(1:EEG.nbchan));
nEegChannels    = sum(eegChannelInd);
x               = EEG.data(eegChannelInd,:)';
disp([num2str(nEegChannels),' EEG channels found']);
EEG             = pop_select(EEG,'channel',eegChannelPos);
Fe              = EEG.srate;

%% Cut pattern
nChan           = EEG.nbchan;
nChanPattern    = round(patternChanProp/100*nChan);
chanPos         = randperm(nChan);
patternChanPos  = chanPos(1:nChanPattern);

patternTimePos  = (1+patternStart*Fe):((1+patternStart*Fe)+patternDuration*Fe);
patternData     = EEG.data(patternChanPos,patternTimePos);


%% Create stim
nNoiseParts         = sum(nRepets)+1;
maxLengthSec        = sum(patternDuration.*nRepets)+(nNoiseParts*noiseDurationInt(2));
stim                = zeros(nChan,maxLengthSec*Fe);
patternInd          = zeros(maxLengthSec*Fe,1);
segPos              = 1;
noisePartStartPos   = round(1+(patternStart+patternDuration)*Fe);
for i=1:nRepets+nNoiseParts
    %- Alternate between noise (not repeated) and patterns
    %- "Noise"
    if rem(i,2)==1
        noiseDuration       = noiseDurationInt(1)+noiseDurationInt(2)*rand(1);
        segmentTimePos      = noisePartStartPos:fix((noisePartStartPos+noiseDuration*Fe));
        segment             = EEG.data(:,segmentTimePos);
        segmentPattern      = 0*ones(1,size(segment,2));
        noisePartStartPos   = noisePartStartPos+length(segmentTimePos);
    %- Pattern (plus noise)
    else
        segmentTimePos      = noisePartStartPos:fix((noisePartStartPos+patternDuration*Fe));
        segment             = EEG.data(:,segmentTimePos);
        %- replace appropriate channels with the pattern
        segment(patternChanPos,:)   = patternData;
        segmentPattern              = 1*ones(1,size(segment,2));
        noisePartStartPos           = noisePartStartPos+length(segmentTimePos);
    end
    stim(:,segPos:segPos+size(segment,2)-1)     = segment;
    patternInd(segPos:segPos+size(segment,2)-1) = segmentPattern;
    segPos = segPos+length(segment);
end
stim        = stim(:,1:segPos-1);
patternInd  = patternInd(1:segPos-1);


%% Write stim as .edf file and .mat file (and save patternInd)
% eegfilename = ['eegstim_',num2str(nPatterns),'patterns_',num2str(patternChanProp),'percentChan_',num2str(nRepets),'_',num2str(patternDuration*1000),'ms.edf'];
% EEGstim             = eeg_emptyset();
% EEGstim.data        = stim;
% EEGstim.chanlocs    = EEG.chanlocs;
% EEGstim.srate       = EEG.srate;
% EEGstim.setname     = eegfilename;
% EEGstim             = eeg_checkset(EEGstim);
% pop_writeeeg(EEGstim,eegfilename,'TYPE','EDF');
% 
% x                   = stim;
% save([eegfilename(1:end-4),'.mat'],'x');
% 
% %- pattern ind
% patternFilename = [eegfilename(1:end-4),'_patternInd.mat'];
% save(patternFilename,'patternInd');

%% Visualize stim
figure;
ax(1) = subplot(4,1,1:3);
imagesc(stim); 
axis tight; ylabel('Electrode');
title(['EEG Stimulus - ',num2str(1000*patternDuration),'ms pattern repeated ',num2str(nRepets),' times']);
ax(2) = subplot(4,1,4);
plot(patternInd); 
axis tight; xlabel('time'); title('Pattern');
linkaxes(ax,'x');

%% Plot 1 channel
nPatterns       = sum(diff(patternInd)==-1);
patternStarts   = find(diff(patternInd)==1);
patternEnds     = find(diff(patternInd)==-1)-1;

% if min([size(x)])==1; x=x'; end;
iChan   = 1;
t       = linspace(0,size(stim,2)/Fe,size(stim,2));
figure; hold on;
plot(t,stim(iChan,:))
for i=1:nPatterns
    patternIndi = patternStarts(i):patternEnds(i);
	plot(t(patternIndi),stim(iChan,patternIndi),'r');
end
axis tight;
xlabel('time (s)'); ylabel('Amp (uV)'); title(['EEG stim - channel ',EEG.chanlocs(iChan).labels]);


