function [spkTimeInd, spkChanInd, vNeuron, Ycomp, Yfilter] = ...
    signal2spikeconverter (x,Fe,nFilters,filterOrder,fMinHz,fMaxHz, ...
    freqScale,rectifierType,GainiSyn,filterOrderLP,freqCharRatio)
% [spkTimeInd, spkChanInd, vNeuron, Ycomp, Yfilter] = 
%    SIGNAL2SPIKECONVERTER (x, Fe, nFilters,filterOrder, fMinHz, fMaxHz,
%    freqScale, rectifierType, GainiSyn, filterOrderLP, freqCharRatio)
%
% Converts the input signal x to a spikes array using a cochlea model.
%
% Pipeline : 
%   1. Normalization of the signal between -1 and 1
%   2. Filter design using the specified frequency scale freqScale
%   3. Filter the input signal x through the filterbank and rectify it 
%   (mexfile: cochlea_bankfilterrectifier_***) given the specified
%   rectifier type
%   4. Compress the range of the signal : Ycomp = Yrect.^(1/3)
%   5. Leaky Integrate and Fire (LIF) model. Input synaptic current is
%   computed by multiplying Ycomp by a the gain GainiSyn
%   (mexfile: cochlea_LIFmodel)
%
% INPUTS:
%   - x             : Input signal
%   - Fe            : Sampling frequency of input signal
%   - nFilters      : Number of filters in the cochlea filterbank
%   - filterOrder   : Filter order (must be a multiple of 2)
%   - fMinHz        : Minimal frequency of analysis in the cochlea
%   - fMaxHz        : Maximal frequency of analysis in the cochlea
%   - freqScale     : Frequency scale for filterbank. Can be 'erb',
%                     'linear', 'eeg' or 'ieeqEqualPower'
%   - rectifierType : Rectifier type. Can be 'classic', 'classif_v2' or
%                     'hybrid'
%   - GainiSyn      : Synaptic gain. Output of compression step is
%                     multiplied by this gain before feeding the LIF neurons
%   - filterOrderLP : Order of the low pass filters used with the
%                     'classic_v2' rectifier
%   - freqCharRatio : Used for 'classic_v2' rectifier
%
% OUTPUTS:
%   - spkTimeInd    : Array of spike time index 
%   - spkChanInd    : Array of spike channel index
%   - vNeuron       : Neurons potential
%   - Ycomp         : Compressed signal
%   - Yfilter       : Filtered signal
%
% Author(s) : Martin Deudon (2016)

%% Normalization of the signal between -1 and 1
a 	= 2./(max(x)-min(x));
b 	= 1-a.*max(x);
x  	= x*diag(a)+repmat(b,size(x,1),1);
disp('Normalizing to [-1,1]');
nChannels = size(x,2);

%% Filter Design
if strcmpi(freqScale,'erb')
    [b,a,freqChar]  = cochlea_erbfilterbank(Fe,fMinHz,fMaxHz,nFilters,filterOrder,0);
elseif strcmpi(freqScale,'linear')
    [b,a,freqChar]  = cochlea_linearfilterbank(Fe,fMinHz,fMaxHz,nFilters,filterOrder,0);
elseif strcmpi(freqScale,'eeg')
    [b,a,freqChar]  = cochlea_eegfilterbank(Fe,filterOrder,1);
    nFilters        = length(freqChar);
elseif strcmpi(freqScale,'ieegEqualPower')
    [b,a,freqChar]  = cochlea_ieegequalpowerfilterbank(Fe,fMinHz,fMaxHz,nFilters,filterOrder,0);
else
    disp('You must specify freqScale (erb or linear)');
end

Ycomp   = zeros(size(x,1),nChannels*nFilters);
Yfilter = zeros(size(x,1),nChannels*nFilters);
wb_h = waitbar(0,'Filtering, rectification and compression');
for i=1:nChannels
    %% Bank filtering with rectifier
    if strcmp(rectifierType,'classic')
        [bLPfilt,aLPfilt]       = butter(1,2*65/Fe,'low');
        [Yfilteri,Yrect]    	= cochlea_bankfilterrectifier_classic(b,a,bLPfilt,aLPfilt,x(:,i));
    elseif strcmp(rectifierType,'classic_v2')
        aLPfilt = zeros(filterOrderLP+1,nFilters);
        bLPfilt = zeros(filterOrderLP+1,nFilters);
        for iFilt=1:nFilters
            [bLPfilt(:,iFilt), aLPfilt(:,iFilt)] = butter(filterOrderLP,(2/Fe)*freqChar(iFilt)/freqCharRatio);
        end      
        [Yfilteri,Yrect]       = cochlea_bankfilterrectifier_classic_v2(b,a,bLPfilt,aLPfilt,x(:,i));
    elseif strcmp(rectifierType,'hybrid')
        %- Characteristic Period (for the rectifier)
        periodChar     	= round(Fe./freqChar);
        periodChar      = periodChar(:);
        [Yfilteri,Yrect]= cochlea_bankfilterrectifier_hybrid(b,a,x(:,i),periodChar);
    else
        error('You must specify rectifierType (erb or linear)');
    end
    %- Compression
    Yrect(Yrect<0) = 0;
    Yfilter(:,1+(i-1)*nFilters:i*nFilters)  = Yfilteri;
    Ycomp(:,1+(i-1)*nFilters:i*nFilters)  	= Yrect.^(1/3);
    waitbar(i/nChannels,wb_h);
end
close(wb_h);

%% LIF neuron model
disp('LIF model');
iSyn  = GainiSyn*Ycomp;
[vNeuron, spkTimeInd, spkChanInd] = cochlea_LIFmodel(iSyn, Fe);
disp([num2str(length(spkTimeInd)),' spikes']);

end

