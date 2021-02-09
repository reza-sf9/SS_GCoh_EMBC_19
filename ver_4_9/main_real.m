clc
clear
% clearvars -except ssv count
close all

sv_init = .001;       % variance


% load inilial points
fr_vec = [12  25];
sr_vec = [256 512 1024]; 
fr = fr_vec(2); % fre
seg = sr_vec(1); % seg

ch_vec = [4 10 20 32 50 60 64];
ch_num = ch_vec(4); % number of channels
% down sample Y_k
down_r_vec = [1 4 8 16 32];
down_rate = down_r_vec(4);

bw_vec = [2 2.5 4.5 8.5];
bw = bw_vec(2);

str_name = sprintf('cfg_Init_fr%d_seg%d_bw%.1f.mat',fr,seg,bw);
load(str_name)
Yk = cfg_Init.Y_k;


Yk = down_Y_k(Yk, down_rate);


% ch_rand = ceil(tot_ch*rand(1, ch_num));
ch_rand = [1 2 3 4 5 6 7 11 12 13 14 17 18 19 24 25 29 33 34 35 36 37 38 40 41 42 43 48 49 54 55 60];
% ch_rand = (1:ch_num);

y_end = length(Yk);

y_l_vec = [0 .22 .4 .5 .75];
y_u_vec = [1 .2 .45 .65 .9]; 
y_l = 1;
y_u = 1;

coef_y_l = y_l_vec(y_l);
coef_y_u = y_u_vec(y_u);

y_interval = [floor(coef_y_l*y_end)+1 floor(coef_y_u*y_end)];
% y_interval = [4000 4100];


% call fucntion for calculate initial parameters
b_coef =2;
% method_ab = 'a';
method_ab = 'new';
[cfg_Init] = initial_param(Yk, y_interval, ch_rand, method_ab, b_coef);

Y_k = cfg_Init.Y_k;
a_init = cfg_Init.a;
b_init = cfg_Init.b;
L_init = cfg_Init.L;
D_init = cfg_Init.D;
mu_init = cfg_Init.mu;
x0 = cfg_Init.x0;

K = length(Y_k);         % number of all stepes

init_Param.L  = L_init;

init_Param.sv  = sv_init;
init_Param.s_ab= [1 0];
init_Param.mu = mu_init;
init_Param.ab = [a_init b_init];

%% GC (calculate GC value with classical approach [PNAS])

% reza 
% win= ch_num;        % length of window that we want to calculate g_k in there
win= 120;        % length of window that we want to calculate g_k in there
over_lap = floor(ch_num*.8);  % length of overlap window

gc_pnas = gc_gk_pnas_calculator(Y_k, win, over_lap);
h= figure;
plot(gc_pnas);
ylim([0 1]),xlim([1 length(gc_pnas)]);
str_tit_pnas = sprintf('fr = %d - Sample rate = %d - ch = %d - down sample = %d', fr, seg, ch_num, down_rate);
title({'GC measure - pnas appraoch', str_tit_pnas})

% str_save_pnas = sprintf('fr%d_sr%d_ch%d_downR%d.png', fr, sr, ch_num, down_rate);
% saveas(h, str_save_pnas)
%%

ParamaUpdate.sv_update = 0; % 1 update, 0 don't update
ParamaUpdate.mu_update = 0;
ParamaUpdate.ab_update = 1;
ParamaUpdate.L_update  = 0;


Iter = 5;
PARAM= [];
PARAM{1}=init_Param;
sv_update = [];

% inter val of x
step_x = 0.01./2;
x  = (-2: step_x :2).';
y = normpdf(x, x0, sqrt(sv_init));
figure, plot(x, y,'o')


step_interval = 1:length(Y_k);
break_con = 0;
thr_ab = 10^-5;
for iter=1:Iter
    iter
    %% filter-smoother
    [Bayes, pmf_gc, s] = gc_filter_smoother(Y_k, x, PARAM{iter});
    BAYES{iter} = Bayes.smoother;
    PMF_GC{iter} = pmf_gc;
    
    figure, imagesc(step_interval,x,(BAYES{iter}))
    
    %% Parameter Update
    updated_param = gc_parameter_update(PARAM{iter}, ParamaUpdate, Bayes);
    PARAM{iter+1}=updated_param;
    
    tempAB = PARAM{iter+1}.ab;
    A{iter} = tempAB(:, 1);
    B{iter} = tempAB(:, 2);
    
    % break part
    if iter>1
        ab_curr = PARAM{iter+1}.ab;
        a_curr = ab_curr(:, 1);
        b_curr = ab_curr(:, 2);
        
        ab_pre = PARAM{iter}.ab;
        a_pre = ab_pre(:, 1);
        b_pre = ab_pre(:, 2);
        
        % mmse
        e_ab(iter) = immse(a_curr, a_pre) + immse(b_curr, b_pre)
        if e_ab<thr_ab
            break_con = 1;
        end
    end
    
    if break_con==1
        break;
    end
    
    
    %% plot gc 
%     h1 = figure;
%     smoother = Bayes.smoother;
%     imagesc(step_interval, x, smoother);
%     title(['posterior estimated of smoother  - iter = ', num2str(iter)]);
%     xlabel('step')
%     ylabel('x_k')
%     colorbar
%     colormap jet
%     
%     str_save1 = sprintf('xk_%d.jpg', iter);
%     saveas(h1, str_save1);
%     
%     h2 =figure;
%     imagesc(step_interval, s, pmf_gc);
%     title(['estimated value for pmf of GC - iter = ', num2str(iter)]);
%     xlabel('step')
%     ylabel('GC')
%     colorbar
%     colormap jet
%     str_save2 = sprintf('gc_%d.jpg', iter);
%     saveas(h2, str_save2);
    
end

data.xk = BAYES;
data.step_interval = step_interval;
data.x = x;
data.gc = PMF_GC;
data.s = s;
data.gc_pnas = gc_pnas;
data.L = L_init;
data.Y_k = Y_k;

% %% a,b
a_up = zeros(ch_num, iter);
b_up = zeros(ch_num, iter);
for i=1:iter
    ab = PARAM{i}.ab;
    a_up(:,i) = ab(:, 1);
    b_up(:,i) = ab(:, 2);
end

data.a = a_up;
data.b = b_up;

data.A = A;
data.B = B;

data.coef_y_l = coef_y_l;
data.coef_y_u = coef_y_u;

str_save = sprintf('data_paper_sv%.4f_stepX%.4f_b%d_ch%d_fr%d_sr%d_downR%d_Yk_%.2f_%.2f_iter%d.mat',...
    sv_init, step_x, b_coef, ch_num,fr,seg,down_rate,coef_y_l,coef_y_u,iter);
save(str_save, 'data')


