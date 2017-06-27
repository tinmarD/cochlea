clear all;

%% Parameters
% wavFilename     = 'C:\Users\deudon\Desktop\M4\sounds\music\batterie_2.wav';
% wavFilename     = 'C:\Users\deudon\Desktop\M4\_Data\sounds\Gaussian_stim_noisy\stim_SNR10db_1patterns_20_200ms.wav';
wavFilename     = 'C:\Users\deudon\Desktop\M4\_Data\sounds\yes_no.wav';
fMinHz          = 20;
fMaxHz          = 15000;
nFilters        = 1000;
filterOrder     = 4;
%- rectifier type
rectifierType   = 'classic';    % 'classic' or 'hybrid'
% rectifierType   = 'hybrid';     % 'classic' or 'hybrid'
%- filter type
freqScale       = 'erb'; % 'linear' or 'erb'
saveResults     = 1;
resultsFigDir   = 'C:\Users\deudon\Desktop\M4\_Results\Speech\erb';
spikeSaveDir    = 'C:\Users\deudon\Desktop\M4\_Results\Speech\erb';
GainiSyn        = 120;
%- STDP
N               = 64;   % Number of spikes in a packet
W               = 32;   % Number of weights
displayThreshold= 15;


%% Load signal
[x,Fe]          = audioread(wavFilename,'double');
if size(x,2)~=1
    warning('Work only with mono wave files');
    x = x(:,1);
end

%% Signal to spike conversion
tic;
[spkTimeInd, spkChanInd] = signal2spikeconverter(x,Fe,nFilters,filterOrder,fMinHz,fMaxHz,freqScale,rectifierType,GainiSyn);
toc;

%% Save spike list
if exist([wavFilename(1:end-4),'.mat'])
    %- variable is called patternInd
    load([wavFilename(1:end-4),'.mat']);
    patternInd                  = patternInd(spkTimeInd+1);
    patternInd(patternInd==0)   = NaN;
    spike_list       = [double(spkTimeInd+1),double(spkChanInd+1),patternInd];
else
    spike_list       = [double(spkTimeInd+1),double(spkChanInd+1),nan(size(spkTimeInd,1),1)];
end
tempSep         = regexp(wavFilename,filesep);
filename        = wavFilename(tempSep(end)+1:end-4);
mkdir(fullfile(spikeSaveDir,filename));
spikelistPath   = fullfile(spikeSaveDir,filename,[filename,'_spike_list_0_001.mat']);
save (spikelistPath,'spike_list');

%% STDP 
tic;
output_spike_list = STDP_spike_packet_v3(spike_list,N,W,nFilters);
toc;

%%
STDP_plot_results_v3 (spike_list, output_spike_list, Fe, nFilters, nFilters, displayThreshold, length(x), 0.010);


%% Results
figure;
plot(spkTimeInd/Fe, spkChanInd,'.','Color',[0.08,0.17,0.55],'MarkerSize',4);
xlabel('Time (s)'); ylabel('Channel'); title('Channel-temporal spike pattern');
axis tight;
