function [A0, An, Bn]= FTtrig_MITlaos(f)
% ============================================
%
% Version 2.1
% Authors: R. H. Ewoldt and P.B. Winter 
% Contact: MITlaos@mit.edu
% Date: 02-Jul-2007
%
% (c) 2007
%
% Please Do Not Distribute,
% contact MITlaos@mit.edu to request
%
% About File:
% Find trigonometric Fourier Series components from FFT:
% f = A0 + SUM_n( An*cos(n*2*pi*t/T + Bn*sin(n*2*pi*t/T)
%
%
% [A0, An, Bn]= FTtrig(f)
%
% VARIABLES
%   f           vector to be transformed
%   A0          essentially mean(f)
%   An          cosine terms
%   Bn          sine terms
%
% SEQUENCE
%   force input to have EVEN number of data points (reqd for fft.m)
%   take FFT > complex vector results
%   extract trigonometric terms from complex vector
%
%
% ============================================


if int16(length(f)/2) == length(f)/2
    %do nothing
else
    %trim last data point to force even number of data points
    %f MUST HAVE EVEN NUMBER OF DATA POINTS!
    d=f; % d is placeholder
    clear f
    f = d(1:length(d)-1);
    clear d
end

n=length(f);  
N=n/2;       %N will be the number of harmonics to consider

Fn(n)=zeros; %initialize complex transform vector
             % this will make it a ROW vector
             % which is necessary for combination later

sizef = size(f);
if sizef(1) == 1 %if ROW vector
    Fn = fft(f); %let Fn be ROW vector
else
    frow = f';
    Fn = fft(frow); %compute FFT Fn=[ low > high | high < low ]
                     %force input to fft to be ROW fector
end
Fn_new = [conj(Fn(N+1))  Fn(N+2:end)  Fn(1:N+1)]; 
%rearrange values such that: Fn_new = [ high < low | low > high ]

Fn_new(:) = Fn_new(:)/n;  %scale results
    
A0 = Fn_new(N+1);
An = 2*real(Fn_new(N+2:end));    %cosine terms
Bn = -2*imag(Fn_new(N+2:end));   %sine terms

