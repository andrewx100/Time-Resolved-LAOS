function [istrain, istress, N, istart, istop]= cycletrim_MITlaos(gamma,tau)
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
% Extract integer numer of strain and stress data, starting with sine
% strain if possible
%
%
% [istrain, istress, N]= cycletrim(gamma,tau)
%
% VARIABLES
%   gamma       strain data input, arbitrary length
%   tau         stress data input, arbitrary length
%   istrain     (integer # of cycles)
%   istress     (aligned with strain)
%   N           # of cycles contained
%
% SEQUENCE
%   (fit sine wave to strain data) - not yet implemented
%   find zeros in strain signal
%   identify first starting point of sine wave
%   identify final point of sine wave
%   output trimmed data
% ============================================

%%% Diagnostic testing initialization
%{
clear 
clc
directory = 'C:\Research\ARES\Janmey\';
files = dir(strcat(directory,'*.txt')) ; 
data=load(strcat(directory,files(8).name)); 
norm = data(500:end,1);  %Normal force data
gamma = data(500:end,2)/100;  %Strain data is expected to be units of percent, so force it to be unitless
tau = data(500:end,3);    %Stress data
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% % % d_zero=[0,0];    %%This was initialized this way to be used with the
% check below, but this has now been superceded by a different process
d_zero=[];

k=0;  %k is a counter for the number of times gamma changes sign
sign_gam = sign(gamma);
for i = 1 : length(gamma)-1
    if sign_gam(i) ~= sign_gam(i+1)
        k=k+1;
        d_zero(k)=i+1;  %index location after sign change
    end
end


%%%%%%% This check is now superceded by the process below, RHE 06.25.2007
% % % if d_zero(1)==0 || d_zero(2)==0
% % %     % if the range selected is too small for there to be less than two sign
% % %     % changes. two sign changes are needed to run the rest of cycletrim 
% % %     % since it is impossible to use a empty array as a placeholder for the  
% % %     % value in the gamma array. 
% % %     
% % %     % give an output before exiting
% % %     istrain = 0;
% % %     istress = 0;
% % %     N       = 0;
% % %     return
% % % end




% if the sign changed, and the previous point was NEGATIVE (or ZERO)
% then it's the beginning of the sine wave

% NB: the following assumes that strain signal is smooth enough so that the
% cycles can be trimmed by finding the gamma=0 crossover points.

lgth = length(d_zero);

if lgth <= 1   
    %if there are 0 or 1 locations of gamma crossing zero,
    %it is impossible to extract the minimum of 1 cycle
  
    % give an output before exiting
    istrain = 0;
    istress = 0;
    N       = 0;
    istart  = 0;
    istop   = 0;
    return
end

if lgth == 2
    %if there are 2 locations where gamma crosses zero,
    %a fancy cycle trimming must be performed, which will NOT start with a
    %sine wave
    %SEQUENCE:  estimate points per cycle
    %           check if there are enough points for single cycle (Npts) (error if not)
    %           include final Npts of signal

    Npts = (d_zero(2) - d_zero(1) ) *2; %estimate number of points per cycle
    
    if length(gamma) < Npts % if there are not enough points
        % give an output before exiting
        istrain = 0;
        istress = 0;
        N       = 0;
        istart  = 0;
        istop   = 0;
        return
    else %if there are enough points
        istart = length(gamma) - Npts + 1;
        istop  = length(gamma);
        N = 1;
    end
end

if lgth > 2
    if lgth/2 ~= round(lgth/2)   % Check if lgth is odd
        istart = d_zero(1);
        istop  = d_zero(end) - 1;
        N = (lgth - 1)/2;
    else
        istart = d_zero(1);
        istop  = d_zero(lgth-1) - 1;
        N = (lgth - 2)/2;
    end
end


%%%%%%%%%%%  Following processing is obsolete, RHE 06.25.2007
% % % % % if gamma(d_zero(1)-1) <= 0
% % % % %     istart = d_zero(1);
% % % % %     init = 1;  %init is used to count the total number of cycles
% % % % % else if gamma(d_zero(2)-1) <= 0
% % % % %         istart = d_zero(2);
% % % % %         init = 2;  
% % % % %     else
% % % % %         errordlg('Something is wrong with istart')
% % % % %     end
% % % % % end
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % if gamma(d_zero(lgth)-1) < 0
% % % % %     istop = d_zero(lgth)-1;
% % % % %     fin = lgth;
% % % % % else if gamma(d_zero(lgth-1)-1) < 0
% % % % %         istop = d_zero(lgth-1)-1;
% % % % %         fin = lgth-1;
% % % % %     else
% % % % %         errordlg('Something is wrong with istop')
% % % % %     end
% % % % % end
% % % % % N = (fin-init)/2;

istrain = gamma(istart:istop);
istress = tau(istart:istop);

% saves the istart and istop values to be used in data preprocessing
% LAOSaddress = getappdata(0 , 'LAOSaddress');
% hdataRange  = getappdata(LAOSaddress , 'hdataRange');
% setappdata(hdataRange  , 'istart', istart);
% setappdata(hdataRange  ,  'istop',  istop);
%%% RHE - removed these lines, since they require hdataRange to be calling
%%% this function (but other things might call this function too!)


if N==round(N)
    %Do nothing
else
    errordlg('Something is wrong with N')
end















