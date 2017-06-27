mex -O CFLAGS="\$CFLAGS -std=c99" cochlea_LIFmodel.c


%% Parameters
Fe              = 44100;
xDurationSec    = 0.2;
fSinHz          = 30; 

%%
t               = linspace(0,xDurationSec,xDurationSec*Fe);
%- Echelon
iSynEchelon     = zeros(xDurationSec*Fe,1);
iSynEchelon(1+round(0.5*xDurationSec*Fe):end) = 1;
%- Click de 1ms
iSynClick1      = zeros(xDurationSec*Fe,1);
iSynClick1(round(0.5*xDurationSec*Fe):round(0.501*xDurationSec*Fe)) = 1;
%- Click de 10ms
iSynClick10     = zeros(xDurationSec*Fe,1);
iSynClick10(round(0.5*xDurationSec*Fe):round(0.51*xDurationSec*Fe)) = 1;
%- Click de 10ms
iSynClick100    = zeros(xDurationSec*Fe,1);
iSynClick100(round(0.5*xDurationSec*Fe):round(0.6*xDurationSec*Fe)) = 1;
%- Sinusoide 
iSynSin                 = sin(2*pi*fSinHz*t);
iSynRect                = iSynSin;
iSynRect(iSynRect<0)    = 0;
[bLPfilt,aLPfilt]       = butter(1,2*65/Fe,'low');
iSynRect                = filter(bLPfilt,aLPfilt,iSynRect);


iSyn = 80*iSynEchelon;
iSyn = iSyn(:);
[vNeuron, spkTimeInd, spkChanInd]       = cochlea_LIFmodel(iSyn, Fe);
[vNeuron2, spkTimeInd2, spkChanInd2]    = cochlea_LIFmodel_2(iSyn, Fe);


%% Results

figure;
ax(1) = subplot(3,1,1);
plot(1000*t,iSyn); title('Input synpatic current');
ylim([min(iSyn)-5,max(iSyn)+5]);
xlabel('time (ms)');
ax(2) = subplot(3,1,2:3);
plot(1000*t,vNeuron); hold on;
% plot(1000*t,vNeuron2,'r--')
plot(xlim,[20,20],'g--');
ylim([-33,40]);
xlabel('time (ms)');
title('Neuron potential');
legend({'With refractory period','With negative reset potential','Spiking threshold'});
linkaxes(ax,'x')
% xlim([4,6]);