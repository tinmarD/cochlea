function [ cf ] = erbspace(lowFreq, highFreq, N, Qear, BWmin)
% [ cf ] = ERBSPACE(lowFreq, highFreq, N, Qear, BWmin)
% This function computes an array of N frequencies uniformly spaced between
% lowFreq and highFreq on an ERB scale. 
%
% For a definition of ERB, see Moore, B. C. J., and Glasberg, B. R. (1983).
% "Suggested formulae for calculating auditory-filter bandwidths and
% excitation patterns," J. Acoust. Soc. Am. 74, 750-753.
%
% See also cochlea_erbfilterbank

if nargin==3
    Qear    = 9.26449;
    BWmin   = 24.7;
end

lowFreq     = double(lowFreq);
highFreq    = double(highFreq);

cf      = -(Qear*BWmin) + exp((N-1:-1:0) * (-log(highFreq + Qear * BWmin) + log(lowFreq + Qear * BWmin))/(N-1)) * (highFreq + Qear * BWmin);

end

