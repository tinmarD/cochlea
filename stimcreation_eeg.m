% STIMCREATION_EEG
% Create a artificial EEG signal with repeated pattern(s) from a real EEG
% file. The eeg part selected for the pattern is fixed with the
% patternStart and the patternDuration variables.
% There can be multiple different patterns. The discontinuities between
% pattern parts and "noise" is removed by adding an offset to the entire
% pattern.
% The pattern appear on a fixed percentage of the channels
% (patternChanProp). These channels are selected randomly.
% 
%
% Author(s) : Martin Deudon (2016)

rng('shuffle');

%% Parameters
inputSignalPath     = 'C:\Users\deudon\Desktop\SAB\_Data\COG_027\SAB_Bipolaire\SAB 600 ms\SAB REC 1.edf';
outputDirPath       = 'C:\Users\deudon\Desktop\M4\_Data\EEGstim';
patternStart        = 90;                       % Starting time of the pattern (s)
patternDuration     = 0.05;                     % Patterns duration (s)
noiseDurationInt    = [0.1 0.2];                % Interval between 2 patterns (s)
patternChanProp     = 50;                       % Proportion of channels where the pattern is present (%)
nRepets             = 50;                      % Number of repetition
patternColors       = {'c','g','r','b','y'};    % Color of each different pattern for visualization
nPatterns           = length(patternDuration);  % Number of different patterns

%% Load signal
EEG = pop_biosig(inputSignalPath);
%- Remove non-eeg channels
eegChannelInd   = cellfun(@(x)strcmp(x(1:3),'EEG'),{EEG.chanlocs.labels});
eegChannelPos   = nonzeros(eegChannelInd.*(1:EEG.nbchan));
nEegChannels    = sum(eegChannelInd);
disp([num2str(nEegChannels),' EEG channels found']);
EEG             = pop_select(EEG,'channel',eegChannelPos);
Fe              = EEG.srate;

%% Cut pattern
nChan           = EEG.nbchan;
nChanPattern    = round(patternChanProp/100*nChan);
disp(['Pattern appears on ',num2str(nChanPattern),'/',num2str(nChan),' channels (',num2str(patternChanProp),'%)']);
patternChanPos  = randperm(nChan,nChanPattern);

patternTimePos  = (1+patternStart*Fe):((1+patternStart*Fe)+patternDuration*Fe);
patternData     = EEG.data(patternChanPos,patternTimePos);


%% Create stim
nNoiseParts         = sum(nRepets)+1;
maxLengthSec        = sum(patternDuration.*nRepets)+(nNoiseParts*noiseDurationInt(2));
stim                = zeros(nChan,ceil(maxLengthSec*Fe));
patternInd          = zeros(ceil(maxLengthSec*Fe),1);
segPos              = 1;
noisePartStartPos   = round(1+(patternStart+patternDuration)*Fe);
for i=1:nRepets+nNoiseParts
    %- Alternate between noise (not repeated) and patterns
    %- "Noise"
    if rem(i,2)==1
        noiseDuration       = noiseDurationInt(1)+rand(1)*diff(noiseDurationInt);
        segmentTimePos      = noisePartStartPos:fix((noisePartStartPos+noiseDuration*Fe));
        segment             = EEG.data(:,segmentTimePos);
        if i>1
            % Offset to remove discontinuities
            segment             = segment - repmat(segment(:,1)-stim(:,segPos-1),1,size(segment,2));
        end
        segmentPattern      = 0*ones(1,size(segment,2));
        noisePartStartPos   = noisePartStartPos+length(segmentTimePos);
    %- Pattern (plus noise)
    else
        segmentTimePos      = noisePartStartPos:fix((noisePartStartPos+patternDuration*Fe));
        segment             = EEG.data(:,segmentTimePos);
        %- replace appropriate channels with the pattern
        segment(patternChanPos,:)   = patternData;
        % Offset to remove discontinuities
        segment                     = segment - repmat(segment(:,1)-stim(:,segPos-1),1,size(segment,2));
        segmentPattern              = 1*ones(1,size(segment,2));
        noisePartStartPos           = noisePartStartPos+length(segmentTimePos);
    end
    stim(:,segPos:segPos+size(segment,2)-1)     = segment;
    patternInd(segPos:segPos+size(segment,2)-1) = segmentPattern;
    segPos = segPos+size(segment,2);
end
stim        = stim(:,1:segPos-1);
patternInd  = patternInd(1:segPos-1);


%% Write stim as .edf file and .mat file (and save patternInd)
eegfilename = ['eegstim2_',num2str(nPatterns),'patterns_',num2str(patternChanProp),'percentChan_',num2str(nRepets),'_',num2str(patternDuration*1000),'ms.edf'];
EEGstim             = eeg_emptyset();
EEGstim.data        = stim;
EEGstim.chanlocs    = EEG.chanlocs;
EEGstim.srate       = EEG.srate;
EEGstim.setname     = eegfilename;
EEGstim             = eeg_checkset(EEGstim);
outputFilePath      = createuniquefilepath(fullfile(outputDirPath,eegfilename));
pop_writeeeg(EEGstim,outputFilePath,'TYPE','EDF');

x                   = stim;
save([outputFilePath(1:end-4),'.mat'],'x');

%- pattern ind
patternFilename = [outputFilePath(1:end-4),'_patternInd.mat'];
save(patternFilename,'patternInd');

disp(['Stim saved in ',outputFilePath]);
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




