function [M,L,EtaM,EtaL,NE,NV]=MITLaos(fi,fs,PPQC,gamma,tauxy)

d_zero=[];

k=0;  %k is a counter for the number of times gamma changes sign
sign_gam = sign(gamma);
for i = 1 : length(gamma)-1
    if sign_gam(i) ~= sign_gam(i+1)
        k=k+1;
        d_zero(k)=i+1;  %index location after sign change
    end
end

% if the sign changed, and the previous point was NEGATIVE (or ZERO)
% then it's the beginning of the sine wave

% NB: the following assumes that strain signal is smooth enough so that the
% cycles can be trimmed by finding the gamma=0 crossover points.

lgth = length(d_zero);
if lgth == 0    
    % give an output before exiting
    gam = 0;
    tau = 0;
    Ncycles = 0;
    errordlg('Selected data never crosses zero.  Unable to find integer number of cycles'         ,'Error')
    return
    
elseif lgth == 1   
    %if there are 0 or 1 locations of gamma crossing zero,
    %it is impossible to extract the minimum of 1 cycle
%     msgbox({'It looks like you have only one cycle, but this cannot be confirmed.';'';'Good luck.'} ,'Warning')
    Ncycles = 1;
    istart  = 1;
    istop   = length(gamma);
    
    gam     = gamma(istart:istop);
    tau     = tauxy(istart:istop);
    
%     if timecol ~=0
%         time = time(istart:istop);
%         setappdata(LAOSaddress,'time',time);
%     end
    
elseif lgth == 2
    %if there are 2 locations where gamma crosses zero,
    %a fancy cycle trimming must be performed, which will NOT start with a
    %sine wave
    %SEQUENCE:  estimate points per cycle
    %           check if there are enough points for single cycle (Npts) (error if not)
    %           include final Npts of signal

    Npts = (d_zero(2) - d_zero(1) ) *2; %estimate number of points per cycle
    
    if length(gamma) < (0.95*Npts) % if there are not enough points
        % give an output before exiting
%         errordlg({'It looks like there aren''t enough points for a complete cycle.';'';'Proceed with extreme caution.'},'Warning')
        istart  = 1;
        istop   = length(gamma);
        Ncycles       = 1;
        
        gam     = gamma(istart:istop);
        tau     = tauxy(istart:istop);
        
%         if timecol ~=0
%             time = time(istart:istop);
%             setappdata(LAOSaddress,'time',time);
%         end

        
    elseif length(gamma) > (1.05*Npts) % check for excess data beyond one cycle (x% tolerance)
        istart = length(gamma) - Npts + 1;
        istop  = length(gamma);
        Ncycles = 1;
        
        gam     = gamma(istart:istop);
        tau     = tauxy(istart:istop);
        
%         if timecol ~=0
%             time = time(istart:istop);
%             setappdata(LAOSaddress,'time',time);
%         end

        msgbox('The beginning of your data was trimmed in order to have exactly one cycle' ,'Warning')
    else
        istart = 1;
        istop  = length(gamma);
        Ncycles = 1;
        
        gam     = gamma(istart:istop);
        tau     = tauxy(istart:istop);

%         if timecol ~=0
%             time = time(istart:istop);
%             setappdata(LAOSaddress,  'time'      ,        time);
%         end
        
%         msgbox('It looks like you have exactly one cycle.  It will not be trimmed.' ,'Warning')
    end
    

    
elseif lgth > 2
    %perform cycle trimming as usual
    [gam, tau, Ncycles, istart, istop] = cycletrim_MITlaos(gamma, tauxy);
%     if timecol ~=0
%         time = time(istart:istop);
%         setappdata(LAOSaddress,  'time'      ,        time);
%     end
% 
end
   
[A0,AnS,BnS]=FTtrig_MITlaos(tau);
[gA0,gAnS,gBnS]=FTtrig_MITlaos(gam);

% gam_0 = abs(gBn(Ncycles)); %identify strain amplitude
gam_0 = sqrt( gBnS(Ncycles)^2 + gAnS(Ncycles)^2 );   % acknowledge possible phase shift, but neglect h.o.t.
delta = atan( gAnS(Ncycles) / gBnS(Ncycles))    ;    % raw signal phase shift


for q=1:length(AnS)  %Create NOT SHIFTED Fourier coefficients
    An(q) = AnS(q)*cos(q*delta/Ncycles) - BnS(q)*sin(q*delta/Ncycles);
    Bn(q) = BnS(q)*cos(q*delta/Ncycles) + AnS(q)*sin(q*delta/Ncycles);
end
% Now proceed as it was before (Introduction of UNshift in v3p3 RHE)  %
% save('D:\Dropbox\Data\Rheology\ARES\2011\06-29-2011\An2.txt','An','-ascii');
% save('D:\Dropbox\Data\Rheology\ARES\2011\06-29-2011\Bn2.txt','Bn','-ascii');


%%% RHE, 06.25.2007, updated now that cycletrim can output a trimmed strain
%%% signal which does not necessarily begin as a positive sine wave.
% Ensure that fundamentals of An and Bn are postive (only check the largest
% coefficient, in case the material is predominantly elastic or viscous, in
% which case one of the coefficient should be very close to zero and is
% therefore not a reliable check
if abs(An(Ncycles)) > abs(Bn(Ncycles))
    if An(Ncycles) < 0
        An = -An;
        Bn = -Bn;
    end
else % if the fundamental of Bn was larger, use it as the reference
    if Bn(Ncycles) < 0
        An = -An;
        Bn = -Bn;
    end
end
%%% RHE, 06.25.2007, updated now that cycletrim can output a trimmed strain
%%% signal which does not necessarily begin as a positive sine wave.

% Determine the sampling frequency of the reconstructed wave


PPQC = 300;

PPC=4*PPQC; %Points Per Cycle
% TP=6*PPQC+1; %Total Points: Points for Three Half-Cycles plus one for overlap

gam_recon=zeros(PPC,1);
for q=1:PPC %sum for each point in time, 1 full cycle, no overlap
    gam_recon(q) = gam_0*sin(2*pi*(q-1)/PPC); %Reconstruct WITHOUT phase shift
end
%make gam_recon 1.5 cycles with 1 point overlap
gam_recon = [gam_recon; gam_recon(1:2*PPQC+1)];

% % w (omega) is currently a MANUAL input
% % strain-rate is equal to omega*strain-shifted-1/4-cycle
% gamdot_recon=w*gam_recon(PPQC+1:PPQC+PPC);  %One cylce of gamdot
% gamdot_recon=[gamdot_recon; gamdot_recon(1:2*PPQC+1)]; %make 1.5 cycles
% % figure;plot(gamdot_recon);hold on;plot(gam_recon)
% 
% %%% Obsolete, now that I'm reconstructing strain data too.
% % T=round(length(gam)/Ncycles) %integer number of data points per cycle
% % Tnum  = length(gam)/Ncycles; %decimal number of data points per cycle

tau_recon = zeros(PPC,1); %initialize tau_recon   (m harmonics included)
% FTtau_e = zeros(PPC,1);
% FTtau_v = zeros(PPC,1);
% 
% tau_recon1 = zeros(PPC,1); %initialize tau_recon1 (1st harmonic only)
% tau_recon3 = zeros(PPC,1); %initialize tau_recon3 (1st & 3rd Harmonics)
% tau_e_3    = zeros(PPC,1); %"elastic" stress, from 1st & 3rd Harmonics
% tau_v_3    = zeros(PPC,1); %"viscous" stress, from 1st & 3rd Harmonics

max_Nr_harm = int32(fs/(2*fi)); % maximum number of harmonics to be consider
m=double(max_Nr_harm - mod(max_Nr_harm,3)); 


for q=1:PPC %sum for each point in time for 1 full cycle, no overlap
    for p=1:2:m %m:total number of harmonics to consider
                %sum over ODD harmonics only
        tau_recon(q) = tau_recon(q) + Bn(Ncycles*p)*sin(p*2*pi*(q-1)/PPC) ...
            + An(Ncycles*p)*cos(p*2*pi*(q-1)/PPC);
        
%         FTtau_e(q)   = FTtau_e(q) + Bn(Ncycles*p)*sin(p*2*pi*(q-1)/PPC);
%         FTtau_v(q)   = FTtau_v(q) + An(Ncycles*p)*cos(p*2*pi*(q-1)/PPC);
% 
    end
%     for p=1:3 %Now just the first 3 harmonics
%         tau_recon3(q) = tau_recon3(q) + Bn(Ncycles*p)*sin(p*2*pi*(q)/PPC) ...
%             + An(Ncycles*p)*cos(p*2*pi*(q)/PPC);
%     end
% %RHE, added June 15, 2007, trying to use FT to reconstruct "Chebyshev"
% %curves
%     for p=1 %Now just the first harmonic
%         tau_recon1(q) = tau_recon1(q) + Bn(Ncycles*p)*sin(p*2*pi*(q-1)/PPC) ...
%             + An(Ncycles*p)*cos(p*2*pi*(q-1)/PPC);
%     end
%     for p=1:2:3 %Harmonics 1 & 3, for "elastic" stress
%         tau_e_3(q) = tau_e_3(q) + Bn(Ncycles*p)*sin(p*2*pi*(q-1)/PPC);
%     end
%     for p=1:2:3 %Harmonics 1 & 3, for "elastic" stress
%         tau_v_3(q) = tau_v_3(q) + An(Ncycles*p)*cos(p*2*pi*(q-1)/PPC);
%     end
end




% %make FTtau_* have one point overlap
% FTtau_e(PPC+1,1)=FTtau_e(1);
% FTtau_v(PPC+1,1)=FTtau_v(1);
% tau_e_3(PPC+1,1)=tau_e_3(1);
% tau_v_3(PPC+1,1)=tau_v_3(1);
% %make tau_recon* 1.5 cycles with 1 point overlap
tau_recon = [tau_recon; tau_recon(1:2*PPQC+1)];
% tau_recon3 = [tau_recon3; tau_recon3(1:2*PPQC+1)];
% figure;plot(gam_recon, tau_recon)

tau_recon_Ncycles = [];
for r=1:Ncycles
    tau_recon_Ncycles = [tau_recon_Ncycles; tau_recon(1:PPC)];
end

% plot(gam_recon,'o')

if Bn(Ncycles)>0  %ensure that G_1' is positive
    Gp = Bn/gam_0;  %G' from sine terms
else
    Gp = -Bn/gam_0;
end

if An(Ncycles) > 0  %ensure that G_1'' is positive
    Gpp = An/gam_0; %G'' from cosine terms
else 
    Gpp = -An/gam_0;
end
M=0;
Lo=0;
EtaM = 0;
EtaL = 0;
for p=1:2:m
    M = M + p*Gp(Ncycles*p);
    Lo = Lo + Gp(Ncycles*p)*(-1)^((p-1)/2);
    
    EtaM = EtaM + (1/2*pi*fi)*p*Gpp(Ncycles*p)*(-1)^((p-1)/2);
    EtaL = EtaL + (1/2*pi*fi)*Gpp(Ncycles*p);
end
%M
L=Lo;

% S=L/M;
% EtaT=EtaL/EtaM;

% S2=(L-M)/L;
NE=(L-M)/Gp(Ncycles);
NV=2*pi*fi*(EtaL-EtaM)/Gpp(Ncycles);
