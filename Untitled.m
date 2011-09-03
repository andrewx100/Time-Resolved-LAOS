filename='D:\Dropbox\Data\Rheology\ARES\2011\06-29-2011\DTimeSwp (4) 000 .txt';
A=importdata(filename);
x=A.data;
[row,col]=size(x);

% Signal and sampling parameters
fs_prime=50;
fi=1;


% Add a time column
t=0:1/fs_prime:(row-1)/fs_prime;
t=t';
x=[t x];

% Total number of cycles
Nr_of_Cycle=int32(t(row)/fi);

% Rows per cycle
row_per_cycle = round(fs_prime/fi);
data=[];

% data_1=x(1:3*row_per_cycle,:);
% [row_1,col_1]=size(data_1);
% 
% freq_bin=0:row_1-1;
% freq_res=fs_prime/row_1;
% freq_1=freq_bin * freq_res;
% Xraw=fft(data_1(:,3))/row_1;
% cutoff=ceil(row_1/2);
% freq_1=freq_1(1:cutoff);
% Xraw=Xraw(1:cutoff);
% plot(freq_1,abs(Xraw))
%     f1 = round(fi / freq_res);
%     f3 = round(3 * fi / freq_res);
%     f5 = round(5 * fi / freq_res);
%     f7 = round(7 * fi / freq_res);
%     f9 = round(9 * fi / freq_res);
%         I31 = abs(Xraw(f3+1))/abs(Xraw(f1+1));
%     I51 = abs(Xraw(f5+1))/abs(Xraw(f1+1));
%     I71 = abs(Xraw(f7+1))/abs(Xraw(f1+1));
%     I91 = abs(Xraw(f9+1))/abs(Xraw(f1+1));




for i = 0:1:Nr_of_Cycle-3;
    begin_row = round(i * fs_prime / fi)+1;
    ti=x(begin_row,1);
    data_i = x(begin_row:begin_row + 2 * row_per_cycle,:); % Take exactly 2 cycles
    
    [row_i, col_i]=size(data_i);
    
    freq_bin = 0:row_i-1; % vector of frequency bins
    freq_res = fs_prime/row_i; % frequency resolution
    freq_i = freq_bin * freq_res; % frequency axis
    Xraw = fft(data_i(:,3))/row_i; % normalized fft
    cutoff = ceil(row_i/2); % only use the first half of the FFT
    freq_i=freq_i(1:cutoff);
    Xraw=Xraw(1:cutoff);
    plot(freq_i,abs(Xraw))
    % frequency bins of the harmonics
    f1 = round(fi / freq_res);
    f3 = round(3 * fi / freq_res);
    f5 = round(5 * fi / freq_res);
    f7 = round(7 * fi / freq_res);
    f9 = round(9 * fi / freq_res);
    
    % Relative intensity of harmonics
    I31 = abs(Xraw(f3+1))/abs(Xraw(f1+1));
    I51 = abs(Xraw(f5+1))/abs(Xraw(f1+1));
    I71 = abs(Xraw(f7+1))/abs(Xraw(f1+1));
    I91 = abs(Xraw(f9+1))/abs(Xraw(f1+1));
    if i==0
        data=[ti I31 I51 I71 I91];
    else
        data=[data;ti I31 I51 I71 I91];
    end
    
end
save([filename 'ft.txt'],'data','-ascii');

