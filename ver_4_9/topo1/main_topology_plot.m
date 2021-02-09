clc
clear
close all


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

% load data of eigenvectors for totpology plott
load('Eig_Info_25')

eig_vec = Eig_Info.eig_vec;
eig_val = Eig_Info.eig_val;

% m(2) = number of time windows that is used to calc GC (Eig vector and value)
m = size(eig_vec);

% acquising max value of EIGENVECTORS of all win times for setting it as
% max value of colorbar
max_eVec_1 = zeros(1 , m(2));
min_eVec_1 = zeros(1 , m(2));

%%% chose 3 points for plotting
last_min = 152;
arbit_min_2_plot = [10 70 130];

last_win = m(2);

coef = last_win/last_min;
vec_eig_arbit = round(arbit_min_2_plot*coef);
%%% 

for count = 1 : length(vec_eig_arbit)
        temp_eig_vec = cell2mat(eig_vec(vec_eig_arbit(count)));
    
    % use diag terms of eigenvector matrix as elements represent
    % eigevectors of each channel corresponding to the biggest eigenvalue
    % because the beiggest eigenvalue is the first eigen value, I used
    % first column of eigenvectors
    first_col_eig_vec = abs(temp_eig_vec(:,1));

    % in addition, the second largest corresponding Eigenvectors of all
    % channels can be assessed, for this aim we should use the second
    % column of eigenvector matrix
    second_col_eig_vec = abs(temp_eig_vec(:,2));
    
    max_eVec_1(1 , count) = max(first_col_eig_vec);
    min_eVec_1(1 , count) = min(first_col_eig_vec);
    
    max_eVec_2(1 , count) = max(second_col_eig_vec);
    min_eVec_2(1 , count) = min(second_col_eig_vec);
end

max_tot_1 = max(max_eVec_1);
min_tot_1 = min(min_eVec_1);

max_tot_2 = max(max_eVec_1);
min_tot_2 = min(min_eVec_1);

% limit of colorbar
colorbar_limit_1 = [min_tot_1 max_tot_1];
colorbar_limit_2 = [min_tot_2 max_tot_2];

% plotting topology of EIGENVECTORS for each window time
for count = 1 : length(vec_eig_arbit)
    disp(vec_eig_arbit(count))
    temp_eig_vec = cell2mat(eig_vec(vec_eig_arbit(count)));
    temp_eig_val = cell2mat(eig_val(vec_eig_arbit(count)));
    
    % using first column  of eigenvector matrix as elements represent
    % eigevectors of each channel corresponding to the largest eigenvalue,
    % because the beiggest eigenvalue is the first eigen value, I used
    % first column of eigenvectors
    first_col_eig_vec = abs(temp_eig_vec(:,1));

    % in addition, the second largest corresponding Eigenvectors of all
    % channels can be assessed, for this aim we should use the second
    % column of eigenvector matrix
    second_col_eig_vec = abs(temp_eig_vec(:,2));
    
    
    % removing ref electrods (I think they are channles number 65 and 66)
    color_map = 'jet';
    str_description_1 = sprintf('Largest * min = %d * EigenValue = %.1f', arbit_min_2_plot(count) , temp_eig_val(1,1));
    ft_plot_topo_rs_2(ch_x, ch_y, first_col_eig_vec , count , colorbar_limit_1 , color_map , str_description_1);
    
%     str_description_2 = sprintf('Second Largest * min = %d * EigenValue = %.1f', arbit_min_2_plot(count) , temp_eig_val(2,2));
%     ft_plot_topo_rs_2(ch_x, ch_y, second_col_eig_vec , count , colorbar_limit_2 , color_map , str_description_2);
%     clc
end

