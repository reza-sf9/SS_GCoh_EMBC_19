function [] = vs_calc_largest_lambda(a, b, xgc, x, t_ss_x)

% figure
% figure, imagesc(step_interval, x, xgc)

x_vec = max(xgc);
lambda = zeros(length(a), length(t_ss_x));

for i =1: length(x_vec)
    
    ind = find(x_vec(i)== xgc(:, i));
    x_max(i) = x(ind);
    
    lambda(:, i) = exp(a+b*x_max(i));
    
    ind_lambda_max(i) = find(lambda(:, i) == max(lambda(:, i)));
    
end

figure
plot(t_ss_x, x_max, 'LineWidth', 3)

hold on
plot(t_ss_x, ind_lambda_max, 'LineWidth', 3)
 
xlabel('time (mins)')

ylim([min(x_max) 3])
xlim([t_ss_x(1) t_ss_x(end)])

legend('estimated x', 'index number of largest \lambda')



end