

%% calcul equal area under the curves between fmin and fMax, when the curve is 
% a function of the type x:->K*x^-b

fmin = 3;
fMax = 75;
nFilters = 6;


%- Function parameters
K = 250;
b = 1.5;
f = @(x)K*x.^(-b);

% Area function (integral)
area = @(x)K*(x.^(-b+1))/(-b+1);

Atot = area(fMax)-area(fmin);

% find the values which divide total Area (Atot) in N equal areas
xFun = @(i,N)(((-b+1)/K)*(i*Atot/N+area(fmin))).^(1/(-b+1));

fVect = xFun(0:nFilters,nFilters);

%% Results

x = fmin:0.01:fMax;
figure; hold on;
plot(x,f(x));

stem(fVect,f(fVect),'r');
axis tight;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title(['iEEG PSD Approximation : x->',num2str(K),'*f^-',num2str(b)]);
