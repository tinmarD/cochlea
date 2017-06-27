%% Parameters
Fe              = 44100;
xDurationSec    = 0.2;
freqSin         = 300; 
freqLP          = 300; 

%%
t               = linspace(0,xDurationSec,xDurationSec*Fe);
%- Echelon
inEchelon     = zeros(xDurationSec*Fe,1);
inEchelon(1+round(0.5*xDurationSec*Fe):end) = 1;
inEchelon2    = 1-inEchelon;
%- Click de 1ms
inClick1      = zeros(xDurationSec*Fe,1);
inClick1(round(0.5*xDurationSec*Fe):round(0.501*xDurationSec*Fe)) = 1;
%- Click de 10ms
inClick10     = zeros(xDurationSec*Fe,1);
inClick10(round(0.5*xDurationSec*Fe):round(0.51*xDurationSec*Fe)) = 1;
%- Click de 100ms
inClick100    = zeros(xDurationSec*Fe,1);
inClick100(round(0.5*xDurationSec*Fe):round(0.6*xDurationSec*Fe)) = 1;
%- Sinusoide 
inSin   = sin(2*pi*freqSin*t);

xIn     = inClick100;
xInRect = xIn;
xInRect(xInRect<0) = 0;

%% 
[bLPfilt,aLPfilt]   = butter(2,2*freqLP/Fe,'low');
yOut                = filter(bLPfilt,aLPfilt,xInRect);
freqz(bLPfilt,aLPfilt);

%% Results
tsec2ind  	= @(x)(x-t(1))*Fe+1;

figure;
plot(t,xInRect); hold on;
plot(t,yOut,'r');
legend({'Input','Filtered signal'});
ylim([-0.2,1.2]);
