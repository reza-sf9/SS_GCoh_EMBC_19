function [] = vs_mtp_specto(ch, T, coef_y_l, coef_y_u , label_size, box_active)
% load('ch_1.mat');

%% cacluating multitaper spectrogram
cfg_mtp_spect = struct;

cfg_mtp_spect.T = T;                % time indices
cfg_mtp_spect.ch_des_num = 1;       % desired channel for plotting spectrogram
cfg_mtp_spect.f_des_l = 1;          % (freq(Hz)) the lower range for plotting spectrogram
cfg_mtp_spect.f_des_u = 30;         % (freq(Hz)) the lower range for plotting spectrogram
cfg_mtp_spect.mtp_NW = 3.5;         % NW (halfbandwidth) value for calculating multitaper
cfg_mtp_spect.mtp_win_length = 64*2;  % (sec) win length for calculating multitaper
cfg_mtp_spect.mtp_over_lap = .75;    % (percentage) length of overlap window based win length


% this fuction plot spectrogram based multitaper approach for desired channel
out = vs_multiTaper_spectrogram(ch , cfg_mtp_spect);

ind_st = round(length(out.t_ind)*coef_y_l); 
if ind_st==0
    ind_st=1;
end
ind_end = round(length(out.t_ind)*coef_y_u); 
XF_DB_des_range_new = out.XF_DB_des_range(ind_st:ind_end, :);
t_ind_new = out.t_ind(1:ind_end-ind_st+1)*60; % convert to mins

%% plot mutlitaper spectrogram
config_subMTP = struct;

config_subMTP.x_start_cordinate = .12;
config_subMTP.x_end_cordinate = .85;
config_subMTP.y_box = .49;
config_subMTP.h_box = .03;
config_subMTP.diff_box_1_2 = .28;
config_subMTP.color_box_1 = [174, 16, 232]./255;
config_subMTP.color_box_2 = [244, 209, 66]./255;
% data for plotting imagesc
config_subMTP.f_ind_des_range = out.f_ind_des_range;
config_subMTP.t_ind_mTaper = t_ind_new;
config_subMTP.XF_DB_des_range = XF_DB_des_range_new;


config_subMTP.xlbl_mtp = 'Freq (Hz)';
config_subMTP.ylbl_mtp = 'Time (mins)';
config_subMTP.lblSize = label_size;

config_subMTP.box_active = box_active;

vs_plot_multitaper2(config_subMTP)

end