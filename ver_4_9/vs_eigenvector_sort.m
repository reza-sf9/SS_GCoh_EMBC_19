function [] = vs_eigenvector_sort(t_ss_x, xk, min_1, min_2, x, a, b, L)

ind_1 = find(t_ss_x>min_1);
ind_1 = ind_1(1);
x_dist_1 = xk(:,ind_1);
ind_x_1 = find(x_dist_1 == max(x_dist_1));
x_max_1 = x(ind_x_1);
D_1 = diag(exp(a+b*x_max_1));
C_1 = L*D_1*L';
[LL_1 , DD_1] = eig(C_1);




ind_2 = find(t_ss_x>min_2);
ind_2 = ind_2(1);
x_dist_2 = xk(:,ind_2);
ind_x_2 = find(x_dist_2 == max(x_dist_2));
x_max_2 = x(ind_x_2);
D_2 = diag(exp(a+b*x_max_2));
C_2 = L*D_2*L';
[LL_2 , DD_2] = eig(C_2);

min_c = min(min(abs(C_1(:)), abs(C_2(:))));
max_c = max(max(abs(C_1(:)), abs(C_2(:))));

h= figure;
subplot(121)
imagesc(abs(C_1))
str_1 = sprintf('min %d            sum of D = %.2f', min_1, sum(D_1(:)));
title(str_1)
colorbar
caxis([min_c max_c])

subplot(122)
imagesc(abs(C_2))
str_2 = sprintf('min %d            sum of D = %.2f', min_2, sum(D_2(:)));
title(str_2)
colorbar
caxis([min_c max_c])

colormap jet

end