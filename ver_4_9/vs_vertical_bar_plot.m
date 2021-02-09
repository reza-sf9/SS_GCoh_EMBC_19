function [] = vs_vertical_bar_plot(vec_initial, vec_final,y_label, label_size)

max_val = max(max(vec_initial), max(vec_final));
min_val = min(min(vec_initial), min(vec_final));

bar_lim = length(vec_final);


% figure('units','normalized','outerposition',[0 0 1 1])
figure
bar(1:bar_lim, [vec_initial vec_final ], 1)
% bar(1:bar_lim, [vec_initial vec_final ], .7, 'r')

xlim([0 bar_lim+1])
set(gca, 'XTick', 1:bar_lim)


ylim([min_val-.5 max_val+.5])
set(gca,'FontSize',label_size)

set(gca, 'XTIck', [8:8:32])


xlabel('Channel','FontSize',label_size)
ylabel(y_label,'FontSize',label_size)

% 
% set(gca, 'YAxisLocation', 'right')

legend('Initial', 'Final')

% set(gcf, 'PaperPosition', [0 0 6 4]);
% print('b_25','-dpng','-r600')

end