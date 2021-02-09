function []= vs_ss_x(time_index, x, xk, label_size, x_interval)

ind_l = find(x>= x_interval(1));
ind_1 = ind_l(1);

ind_2 = find(x>= x_interval(2));
ind_2 = ind_2(1);

x2 = x(ind_1:ind_2);
xk2 = xk(ind_1:ind_2, :);

figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(time_index, x2, xk2);
% title('distribution of x_k');
xlabel('Time (mins)')
ylabel('x')
colorbar
view([0 -90]);

% colormap jet 
c1 = [255, 255, 255]./255;
c2 = [250, 255, 255]./255;
c3 = [246, 225, 255]./255;
c4 = [240, 200, 255]./255;
c5 = [239, 181, 255]./255;
c6 = [233, 150, 255]./255;
c7 = [227, 134, 255]./255;
c8 = [220, 100, 255]./255;
c9 = [217, 83, 255]./255;
c10 = [212, 50, 255]./255;
c11 = [210, 38, 255]./255;
c12 = [205, 20, 250]./255;
c13 = [195, 5, 245]./255;
c14 = [170, 2, 230]./255;
c15 = [168, 0, 212]./255;
c16 = [155, 0, 200]./255;
c17 = [152, 0, 191]./255;
c18 = [100, 0, 150]./255;
c19 = [83, 0, 105]./255;


colormap([c1; c2; c3; c4; c5; c6; c7; c8; c9; c10; c11; c12; c13; c14; c15; c16; c17; c18; c19;]);


set(gca,'FontSize',label_size)

set(gcf, 'PaperPosition', [0 0 6 4]);
% print('X_12','-dpng','-r600')

end