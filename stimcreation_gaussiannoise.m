% STIMCREATION_GAUSSIANNOISE
% Create a gaussian noise stimulus with repeated patterns. There can be
% multiple different patterns. The signal to noise ratio can be adjusted
% with the SNRdb parameter. 
%
% Author(s) : Martin Deudon (2016)

%% Parameters 
resDirPath          = 'C:\Users\deudon\Desktop\M4\_Data\audioStim';
Fe                  = 44100;                    % Sampling frequency (Hz)
patternDuration     = [0.2];                    % Pattern durations (s)
noiseDurationInt    = [0.1 0.5];                % Time interval for the duration between 2 patterns (s)
nRepets             = [20];                     % Number of repetition of each pattern 
nPatterns           = length(patternDuration);  % Number of patterns
patternColors       = {'c','g','r','b','y'};    % Color of each pattern for visualization
SNRdb               = 6;                        % Signal to noise ratio

stdSignal           = 1;
stdNoise            = stdSignal/(10^(SNRdb/20)); % stdSignal must be 1, zero mean
stdSum              = sqrt(stdSignal^2+stdNoise^2);

rng shuffle;    	% Change the seed of the random number generator
%% Create repeated noise patterns
noisePattern = cell(nPatterns);
for i=1:nPatterns
   noisePattern{i} = stdSignal*randn(round(patternDuration(i)*Fe),1);    
end

nNoiseParts     = sum(nRepets)+1;
maxLengthSec    = sum(patternDuration.*nRepets)+(nNoiseParts*noiseDurationInt(2));
minLengthSec    = sum(patternDuration.*nRepets)+(nNoiseParts*noiseDurationInt(1));
stim            = zeros(maxLengthSec*Fe,1);
patternInd      = zeros(maxLengthSec*Fe,1);

%- Create random pattern order
patternOrder    = [];
for i=1:nPatterns
    patternOrder = [patternOrder,i*ones(1,nRepets(i))];
end
patternOrder    = patternOrder(randperm(length(patternOrder)));

segPos = 1;
% Create stim
for i=1:nNoiseParts+sum(nRepets)
    %- Alternate between noise (not repeated) and patterns
    if rem(i,2)==1
        noiseDuration   = noiseDurationInt(1)+noiseDurationInt(2)*rand(1);
        segment         = stdSum*randn(round(Fe*noiseDuration),1);
        segmentPattern  = 0*ones(1,length(segment));
    else
        segment     	= noisePattern{patternOrder(i/2)};
        segment         = segment+stdNoise*randn(length(segment),1); % Add noise to the pattern
        segmentPattern  = patternOrder(i/2)*ones(1,length(segment));
    end
    stim(segPos:segPos+length(segment)-1)       = segment;
    patternInd(segPos:segPos+length(segment)-1) = segmentPattern;
    segPos = segPos+length(segment);
end
stim        = stim(1:segPos);
patternInd  = patternInd(1:segPos);


%% Normalization
a 	= 2./(max(stim)-min(stim));
b 	= 1-a.*max(stim);
stim  	= stim*diag(a)+repmat(b,size(stim,1),1);
disp('Normalizing to [-1,1]');
disp(['Stim standard deviation : ',num2str(std(stim))]);


%% Write file
wavFilename = ['stim_SNR',num2str(SNRdb),'db_',num2str(nPatterns)','patterns_',];
for i=1:nPatterns
    wavFilename = [wavFilename,num2str(nRepets(i)),'_',num2str(patternDuration(i)*1000),'ms_'];
end
wavFilename     = [wavFilename(1:end-1),'.wav'];
wavFilepath     = createuniquefilepath(fullfile(resDirPath,wavFilename));
patternFilepath = [wavFilepath(1:end-4),'.mat'];
audiowrite(wavFilepath,stim,Fe);
% Save pattern id
save(patternFilepath,'patternInd');


%% Visualize stim
t = linspace(0,length(stim)/Fe,length(stim));
figure; 
ax(1) = subplot(211); hold on;
plot(t(patternInd==0),stim(patternInd==0),'b');
for i=1:nPatterns
    plot(t(patternInd==i),stim(patternInd==i),'Color',patternColors{i});
end
xlabel('time (s)');
axis tight;
ax(2) = subplot(212); hold on;
plot(t,patternInd);
axis tight;
plot(xlim,[0,0],'k');
linkaxes(ax,'x');
