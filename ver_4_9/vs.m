clc
clearvars -except data_det T Fs
close all

% % % % % % scalp L
% % % % %
% % % % % vs_topoPlot_sub6(12)
% % % % % vs_topoPlot_sub6(25)

label_size = 20;
%% data
%%% loading data
if exist('data_det') ~=1
    load('D:\Research\ali_new_work\code\Global Coherence\MainCode\eeganes07laplac250_detrend_all.mat');
    
end

% data_12 = load('data_paper_sv0.0010_stepX0.0050_b2_ch50_fr12_sr256_downR16_Yk_0.00_1.00_iter10.mat');
% load('data_paper_sv0.0005_stepX0.0033_b2_ch20_fr25_sr256_downR1_Yk_0.22_0.45_iter10.mat');



load('data_paper_sv0.0010_stepX0.0050_b2_ch32_fr12_sr256_downR16_Yk_0.00_1.00_iter5.mat');
% load('data_paper_sv0.0010_stepX0.0050_b2_ch32_fr25_sr256_downR16_Yk_0.00_1.00_iter5.mat');

xk = data.xk;
step_interval = data.step_interval;
x = data.x;
PMF_GC = data.gc;
gc_pnas = data.gc_pnas;
s = data.s;
a = data.a;
b = data.b;
coef_y_l = data.coef_y_l;
coef_y_u = data.coef_y_u;
A = data.A;
B = data.B;
L = data.L;



% item = 10;
% for i=1:length(B)
%    tempB = B{i};
%    bb(i) = tempB(item);
% end

% figure, plot(bb),title(['b', num2str(item)])
% xlabel('iteration')
% xlim([1 length(bb)])

% time index
% load('ch_1.mat');
T_min = T./(60);



t_end = T_min(end)*(coef_y_u);

if coef_y_l~= 0
  t_start = T_min(end)*coef_y_l;  
else 
  t_start = T_min(1);  
end

step_t_pnas = (t_end - t_start)./length(gc_pnas);
t_pnas_gc = (t_start :step_t_pnas :t_end-step_t_pnas);

step_t_ss_x = (t_end - t_start)./length(step_interval);
t_ss_x = (t_start :step_t_ss_x :t_end-step_t_ss_x);

color_1 = [174, 16, 232]./255;
color_2 = [0, 0, 255]./255;

%% plot largest eigenvcector and how change it during time

larg_ev_plt = 0; 
if larg_ev_plt == 1
    iter = 4;
    vs_calc_largest_lambda(a(:, iter), b(:, iter), xk{iter}, x, t_ss_x)
end

%% empirical GC (pair of frequencies)
emp_gc_calc = 0;

if emp_gc_calc==1
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(t_pnas_gc, data_12.data.gc_pnas, 'LineWidth', 6, 'Color', color_1)
    hold on
    plot(t_pnas_gc, gc_pnas, 'LineWidth', 6, 'Color', color_2)
    xlim([t_pnas_gc(1) t_pnas_gc(end)])
    ylabel('GC', 'FontSize',label_size), xlabel('Time (mins)', 'FontSize',label_size)
    legend('12 Hz', '25 Hz')
    ylim([0 1])
    set(gca,'FontSize',label_size)
end

%% indvidual GC 
emp_gc_indvidual = 0;

if emp_gc_indvidual==1
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(t_pnas_gc, data.gc_pnas, 'LineWidth', 8, 'Color', 'b')
    
    set(gca, 'YTIck', [0:.5:1])
    set(gca, 'YAxisLocation', 'right')
    set(gca, 'XTIck', [20:40:140])
    
    xlim([t_pnas_gc(1) t_pnas_gc(end)])
    ylabel('GCoh', 'FontSize',label_size), xlabel('Time (mins)', 'FontSize',label_size)
    
    ylim([0 1])
    set(gca,'FontSize',label_size)
    
    set(gcf, 'PaperPosition', [0 0 6 4]);
%     print('Gcoh_25','-dpng','-r600')
end


% MTP spectrogram
mtp_calc = 0;

if mtp_calc == 1
    
    box_active = 0; % if it is 1 we have 2 rectangular box if 0 we don't have them 
    
    ch_num1 = 6;
    ch_num2 = 40;
    ch_num3 = 37;
    ch_num4 = 60;
    
    ch_a = data_det(:,ch_num1);
    ch_b = data_det(:,ch_num2);
    ch_c = data_det(:,ch_num3);
    ch_d = data_det(:,ch_num4);
  
   
%     vs_mtp_specto(ch_a, T, coef_y_l, coef_y_u, label_size, box_active)
    vs_mtp_specto(ch_b, T, coef_y_l, coef_y_u, label_size, box_active)
%     vs_mtp_specto(ch_c, T, coef_y_l, coef_y_u, label_size, box_active)
%     vs_mtp_specto(ch_d, T, coef_y_l, coef_y_u, label_size, box_active)
end

%% time plot
time_plt_calc = 0;
if time_plt_calc ==1
    for num_x=1:64
        ch_x = data_det(:,num_x);
        h = figure;
        plot(T_min,ch_x)
        xlim([0 T_min(end)])
        title([num2str(num_x)])
        str_save = sprintf('time_%d.png', num_x);
        saveas(h, str_save);
        close
        %     vs_mtp_specto(ch_x, T, coef_y_l, coef_y_u, label_size, num_x)
    end
end


% pnas gc
color_gc = [69, 69, 216]./255;
% vs_pnas_gc(t_pnas_gc, gc_pnas, color_gc)


% SS gc
ss_gc_plt = 1;

gc_pnas = data.gc_pnas;
t_pnas_gc; 

if ss_gc_plt == 1
    
    % gc_interval = [.2 .8];
    gc_interval = [0 1];
    for i=5:5
        vs_ss_gc(t_ss_x, s, (PMF_GC{i}).^.01, label_size, gc_interval, t_pnas_gc, gc_pnas)
    end
    
end

% SS x
ss_x_plot = 0;

if ss_x_plot == 1
    
    % x_interval = [-.9 .5];
    x_interval = [-2 2];
    for i= 4:4
        vs_ss_x(t_ss_x, x, (xk{i}).^.01, label_size, x_interval)
    end
    
end


% hold on
% plot(t_pnas_gc,gc_pnas,'Color', [1 1 1])

ab_plt = 1;

if ab_plt == 1
    num = 5;
    min_1 = 10;
    min_2 = 20;
    % vs_eigenvector_sort(t_ss_x, xk{num}, min_1, min_2, x, a(:,num), b(:,num), L)
    
    % a plot
    y_label = 'a';
    a_initial = a(:, 1);
    a_final = a(:, 5);
    vs_vertical_bar_plot(a_initial, a_final, y_label, label_size);
    
    % b plot
    y_label = 'b';
    b_initial = b(:, 1);
    b_final = b(:, 5);
    vs_vertical_bar_plot(b_initial, b_final, y_label, label_size);
    
end
