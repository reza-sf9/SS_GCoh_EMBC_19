function []= vs_pnas_gc(t_pnas_gc, gc_pnas, color_gc)

figure;

plot(t_pnas_gc, gc_pnas,'LineWidth', 2, 'Color', color_gc);
ylim([0 1])
xlim([t_pnas_gc(1) t_pnas_gc(end)]);
xlabel('Time (mins)')
ylabel('GC')


end