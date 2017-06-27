function [ fVect ] = iEEGequalpowerbands(fMinHz,fMaxHz,nFilters,showResults)
% [ fVect ] = iEEGequalpowerbands(fMinHz,fMaxHz,nFilters,showResults)
% calcul equal area under the curves between fmin and fMax, when the curve is 
% a function of the type x:->K*x^-b
% Return the limits of the frequency bands of equal power
%
% Author(s) : Martin Deudon (2016)

if nargin==3
    showResults=0;
end

%- Function parameters
K = 250;
b = 1.5;
f = @(x)K*x.^(-b);

% Area function (integral)
area = @(x)K*(x.^(-b+1))/(-b+1);

Atot = area(fMaxHz)-area(fMinHz);

% find the values which divide total Area (Atot) in N equal areas
xFun    = @(i,N)(((-b+1)/K)*(i*Atot/N+area(fMinHz))).^(1/(-b+1));

fVect   = xFun(0:nFilters,nFilters);
fVect   = fVect(:);

%% Results
if showResults
    x = fMinHz:0.01:fMaxHz;
    figure; hold on;
    plot(x,f(x));

    stem(fVect,f(fVect),'r');
%     axis tight;
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title(['iEEG PSD Approximation : x->',num2str(K),'*f^-',num2str(b)]);
end

end

