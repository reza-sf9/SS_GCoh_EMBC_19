function vs_plot_multitaper(config_subMTP,tit)

% Extract data from configuration

box_active = config_subMTP.box_active;

f_ind_des_range = config_subMTP.f_ind_des_range;
t_ind_mTaper = config_subMTP.t_ind_mTaper;
XF_DB_des_range = config_subMTP.XF_DB_des_range;
xlbl_mtp = config_subMTP.xlbl_mtp;
ylbl_mtp = config_subMTP.ylbl_mtp;
lblSize = config_subMTP.lblSize;

x_start_cordinate = config_subMTP.x_start_cordinate;
x_end_cordinate = config_subMTP.x_end_cordinate;
y_box = config_subMTP.y_box;
h_box = config_subMTP.h_box;
diff_box_1_2 = config_subMTP.diff_box_1_2;
color_box_1 = config_subMTP.color_box_1;
color_box_2 = config_subMTP.color_box_2;

h = figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(f_ind_des_range ,t_ind_mTaper , XF_DB_des_range);
% Create xlabel
xlabel(xlbl_mtp,'FontSize',lblSize);
% Create ylabel
ylabel(ylbl_mtp,'FontSize',lblSize);

set(gca,'FontSize',lblSize)

% Uncomment the following line to preserve the X-limits of the axes
xlim([f_ind_des_range(1) f_ind_des_range(end)]);
ylim([t_ind_mTaper(1) t_ind_mTaper(end)]);

view([-90 90]);
% axis('ij');

colormap jet
colorbar

set(gca, 'YTIck', [20:40:140])

set(gcf, 'PaperPosition', [0 0 6 4]);
% print('mtp_40','-dpng','-r600')

% title([num2str(tit)])
% str_save = sprintf('mtp_%d.png', tit);
% saveas(h, str_save);
% close

%% box 1

if box_active == 1

x_box = x_start_cordinate;
w_box = x_end_cordinate - x_start_cordinate;


annotation(h,'rectangle',...
    [x_box y_box w_box h_box],...
    'Color' , color_box_1,...
    'LineWidth' , 4.5);
%% box 2
x_box_2 = x_start_cordinate;
w_box_2 = x_end_cordinate - x_start_cordinate;
h_box_2 = h_box;
y_box_2 = y_box + diff_box_1_2;

annotation(h,'rectangle',...
    [x_box_2 y_box_2 w_box_2 h_box_2],...
    'Color' , color_box_2,...
    'LineWidth' , 4.5);
end

end
