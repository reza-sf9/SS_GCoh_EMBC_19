function out = vs_multiTaper_spectrogram(data_chs , cfg_mtp_spect)
% this function plots spectrogram based multitpaer approach
%
% INPUTS
%  data_chs       : all channels time series information
% cfg_mtp_spect (arrays of this structure are):
%
% T              : time indices
% ch_des_num     : desired channel for plotting multitaper spectrogram 
% f_des_l        : (freq(Hz)) the lower range for plotting spectrogram
% f_des_u        : (freq(Hz)) the upper range for plotting spectrogram
% mtp_NW         : NW (halfbandwidth) value for calculating multitaper
% mtp_win_length : (sec) win length for calculating multitaper
% mtp_over_lap   : (percentage) length of overlap window based win length


T = cfg_mtp_spect.T;              
ch_des_num = cfg_mtp_spect.ch_des_num;     
f_des_l = cfg_mtp_spect.f_des_l;   
f_des_u = cfg_mtp_spect.f_des_u;   
mtp_NW = cfg_mtp_spect.mtp_NW;
mtp_win_length = cfg_mtp_spect.mtp_win_length;
mtp_over_lap = cfg_mtp_spect.mtp_over_lap;




ch_des = data_chs(: , ch_des_num);
m = size(ch_des);
signal_sample = m(1);

NW = mtp_NW;

win_sec = mtp_win_length;
Fs = 256;
win_sample = win_sec*Fs;
overlap_sample = floor(win_sample .* mtp_over_lap);

jump_sample = win_sample - overlap_sample;
 
% getting first and last indices of each window 
vec_step = 1:win_sample;
i = 0;
ex_con = 0;
while (ex_con == 0)
   i = i+1;
    temp_ind = (i-1)*jump_sample + vec_step;
    ind_end = temp_ind(end);
    if ind_end > signal_sample
        ex_con =1;
    end
    temp_ind_overlap(1,i) = temp_ind(1);
    temp_ind_overlap(2,i) = temp_ind(end);
    clear temp_ind
end

ind_overlap = temp_ind_overlap(:,1:end-1);
% ind_overlap has 2 ROWS 
% first row indicates lower index of each overlapped window
% second row indicates upper index of each overlapped window
% COLUMNS shows value of first and last index number for different windows 


% T scale for hours
T = T./(60*60);


TT = zeros(1 , length(ind_overlap));
for i = 1 : length(ind_overlap)
    ind_l = ind_overlap(1,i);
    ind_u = ind_overlap(2,i);
    ind_mid = floor((ind_l + ind_u )./ 2);
    TT(1,i) = T(ind_mid);
end


XF = [];
for i=1 : length(ind_overlap)
    
    disp(i)
    ind_l = ind_overlap(1,i);
    ind_u = ind_overlap(2,i);
    x =  pmtm(ch_des(ind_l:ind_u), NW , Fs);
    XF = [XF;x'];    
end

m = size(XF);
% signal_sample = number of windows
% m(2) = number of freq 
f_ind = 1:m(2);
t_ind = TT;


% MultiTaper Spectrogram matrix
XF_DB = 20*log(XF);
% seprate data of desired fre range
XF_DB_des_range = XF_DB(: , f_des_l : f_des_u);
f_ind_des_range = f_ind(f_des_l : f_des_u);

% Ordinary Spectrogram matrix
[s] = spectrogram(ch_des , win_sample , overlap_sample , Fs);
s_DB = abs(s.');
s_DB = 20*log(s_DB);
% seprate data of desired fre range
S_DB_des_range = s_DB(: , f_des_l : f_des_u);

% range of colorbar
max_mtp_spect = max(XF_DB_des_range(:));
min_mtp_spect = min(XF_DB_des_range(:));

max_ord_spect = max(S_DB_des_range(:));
min_ord_spect = min(S_DB_des_range(:));

max_colorBar = max(max_mtp_spect , max_ord_spect);
min_colorBar = min(min_mtp_spect , min_ord_spect);


out = struct;
out.f_ind_des_range = f_ind_des_range;
out.t_ind = t_ind;
out.XF_DB_des_range = XF_DB_des_range;
out.min_colorBar = min_colorBar;
out.max_colorBar = max_colorBar;
out.ind_overlap = ind_overlap;

end