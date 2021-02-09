clc
% clear
close all

load('data_paper_sv0.010_b2_ch20_fr12_sr256_downR16_iter7.mat')

fr = 12;
sr = 1024;

xk = data.xk;
step_interval = data.step_interval;
x = data.x;
PMF_GC = data.gc;
s = data.s;
a_up = data.a;
b_up = data.b;

ch = length(b_up);

Iter = length(xk);

sub_1 = [1 2 5 6];
sub_2 = [3 4 7 8];
sub_3 = [9 10];
sub_4 = [13 14];
sub_5 = [11 12];
sub_6 = [15 16];

ind_1 = 1;
ind_2 = 3;

for iter=1 : Iter
    close
    h = figure('units','normalized','outerposition',[0 0 1 1]);

    temp_xk = xk{iter};
    temp_gc = PMF_GC{iter};
    
    subplot(4,4, sub_1);
    imagesc(step_interval, x, temp_xk);
    title('distribution of x_k');
%     xlabel('step')
    ylabel('x_k')
%     colorbar 
    colormap jet
    
    subplot(4,4, sub_2);
    imagesc(step_interval, s, temp_gc);
    title('distribution of gc');
%     xlabel('step')
    ylabel('GC')
    colorbar 
    colormap jet
   
    subplot(4,4, sub_3);
    plot(a_up(ind_1, :),'b','LineWidth', 2);
    hold on 
    plot(iter,a_up(ind_1, iter),'o',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b')

    title(['a',num2str(ind_1)]);
    xlim([1 Iter])
    
    subplot(4,4, sub_4);
    plot(a_up(ind_2, :),'b','LineWidth', 2);
    hold on 
    plot(iter,a_up(ind_2, iter),'o',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b')

    title(['a',num2str(ind_2)]);
    xlim([1 Iter])
    xlabel('iteration')
    
    
    subplot(4,4, sub_5);
    plot(b_up(ind_1, :),'r','LineWidth', 2);
    hold on 
    plot(iter,b_up(ind_1, iter),'o',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r')
    
    title(['b',num2str(ind_1)]);
    xlim([1 Iter])
    
%     xlabel('iteration')
    
    subplot(4,4, sub_6);
    plot(b_up(ind_2, :),'r','LineWidth', 2);
    hold on 
    plot(iter,b_up(ind_2, iter),'o',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r')
    
    title(['b',num2str(ind_2)]);
    xlim([1 Iter])
    xlabel('iteration')
    
    
    suptitle(['Iter =',num2str(iter)])
    
    str_save = sprintf('ab_fr%d_sr%d_iter%d_ch%d_WinNum_%d.png',fr,sr,Iter,ch,iter);
    saveas(h, str_save);
end

j