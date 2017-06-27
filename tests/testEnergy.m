
wavFilename     = 'C:\Users\deudon\Desktop\M4\sounds\Gaussian_sounds\RN1.wav';
fMinHz          = 20;
fMaxHz          = 20000;
nFilters        = 300;
filterOrder 	= 4;    % Only multiple of 2 order
freqBandi       = 3;

%% Read file
Fe = 44100;
x = gensound(Fe*1,1,Fe);
% [x,Fe]          = audioread(wavFilename,'double');
if size(x,2)~=1;
    x = x(:,1);
end

%% Normalization of audio file
a       = 2/(max(x)-min(x));
b       = 1-a*max(x);
xNorm   = a*x+b;
x       = xNorm;
disp('Normalizing to [-1,1]');

%% Filter 
[b,a,freqChar,ERB]  = cochlea_erbfilterbank(Fe,fMinHz,fMaxHz,nFilters,filterOrder);
[bLPfilt,aLPfilt]   = butter(1,2*65/Fe,'low');
[Yfilt,Yrect]       = cochlea_bankfilterrectifier_classic(b,a,bLPfilt,aLPfilt,x);

%% Fourier transform
xSel        = Yfilt(:,freqBandi);
Xselfft     = fft(xSel);

%% Energy
Et = sum(abs(xSel).^2);
Ef = sum(abs(Xselfft.^2))/length(x);

EtTot = sum(abs(x).^2);
disp(['Temporal energy for freq. band around ',sprintf('%.2f Hz : ',freqChar(freqBandi)),num2str(Et)]);
disp(['Frequential energy for freq. band around ',sprintf('%.2f Hz : ',freqChar(freqBandi)),num2str(Ef)]);
disp(['Total energy : ', num2str(EtTot)]);
disp(['Total energy normalized per sample: ', num2str(EtTot/length(x))]);

%% Energy evolution
coeff       = 0.5;
yFiltEnergy = sum(abs(Yfilt).^2);
YfiltNorm   = Yfilt./repmat(sqrt(ERB),length(x),1);

figure;
imagesc(Yfilt');
set(gca,'Ydir','normal');
title('Output of bank filtering'); ylabel('Channel'); xlabel('time');
figure;
imagesc(YfiltNorm');
set(gca,'Ydir','normal');
title('Output of bank filtering normalized'); ylabel('Channel'); xlabel('time');
