function [] = vs_topoPlot_sub6(n_target)
clc

load('L_D_64_12hz.mat')
L_12 = L_init;
D_12 = D_init;

load('L_D_64_25hz.mat')
L_25 = L_init;
D_25 = D_init;

L_tot = zeros(64,64,2);
L_tot(:,:,1) = L_12;
L_tot(:,:,2) = L_25;

switch n_target
    case 12
        L_target = L_12;
    case 25
        L_target = 25;
end 
    


cfg = [];

% for differnet layout specification we must change cfg.layout name
% check following website for diffrent layout's name
% http://www.fieldtriptoolbox.org/template/
layout_name = 'easycapM11.mat';
cfg.layout = layout_name;

% with ft_prepare_layout function we can obtain specification of arbitray
% electrode with using layout name and puuting data of that in layout
% folder of fieldtrip toolbox (I downloaded easycapM11.mat from the
% internet-github- and then I put it in layout folder of FIELDTRIP toolbox)
layout_info = ft_prepare_layout(cfg);

layout_rs.fileName = layout_name;
layout_rs.pos = layout_info.pos;
layout_rs.width = layout_info.width;
layout_rs.height = layout_info.height;
layout_rs.label = layout_info.label;

% save structure contains position and label name of layout to use it for
% plotting topology
% % % % % % % str_name = sprintf('layout_%s_rs.mat',layout_name);
% % % % % % % save('layout_easycapM11_rs.mat', 'layout_rs')
% % % % % % % savdir = 'D:\Ali new work\code\Global Coherence\topo_plot';
% % % % % % % save(fullfile(savdir,'layout_easycapM11_rs.mat'),'layout_rs');
% % % % % % %


fileName = layout_rs.fileName;
pos = layout_rs.pos;
label = layout_rs.label;

ch_x = pos(:,1); % x position of channel site in scalp
ch_y = pos(:,2); % y position of channel site on scalp

% removing refrence electrode position
ch_x = ch_x(1:64);
ch_y = ch_y(1:64);





% acquising max value of EIGENVECTORS of all win times for setting it as
% max value of colorbar
max_eVec_1 = zeros(1 , 2);
min_eVec_1 = zeros(1 , 2);
for iter = 1 : 2
    
    L_temp = L_tot(:,:,iter);
    first_col_eig_vec = L_temp(:,1);
    
    second_col_eig_vec = L_temp(:,2);
    
    max_eVec_1(1 , iter) = abs(max(first_col_eig_vec));
    min_eVec_1(1 , iter) = abs(min(first_col_eig_vec));
    
    max_eVec_2(1 , iter) = abs(max(second_col_eig_vec));
    min_eVec_2(1 , iter) = abs(min(second_col_eig_vec));
end
max_tot_1 = max(max_eVec_1);
min_tot_1 = min(min_eVec_1);

max_tot_2 = max(max_eVec_2);
min_tot_2 = min(min_eVec_2);



% limit of colorbar
colorbar_limit_1 = [min_tot_1 max_tot_1];
colorbar_limit_2 = [min_tot_2 max_tot_2];

% plotting topology of EIGENVECTORS for each window time


first_col_eig_vec = abs(L_target(:,1));

second_col_eig_vec = abs(L_target(:,2));


% removing ref electrods (I think they are channles number 65 and 66)
color_map = 'jet';


% structure for plotting topology plot
config_topoPlot = struct;

% these parameter are equal for different plots
config_topoPlot.ch_x = ch_x;
config_topoPlot.ch_y = ch_y;

%% largest
% these parameter are not equal for different plots
config_topoPlot.data = first_col_eig_vec;
config_topoPlot.colorbar_limit = colorbar_limit_1;

figure
vs_ft_plot_topo_rs_2_6(config_topoPlot);

%% second largest

% these parameter are not equal for different plots
config_topoPlot.data = second_col_eig_vec;
config_topoPlot.colorbar_limit = colorbar_limit_2;

figure
vs_ft_plot_topo_rs_2_6(config_topoPlot);

clc
% end

end

