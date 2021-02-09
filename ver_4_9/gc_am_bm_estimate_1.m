function [param_ab_estimate, diff_am, diff_bm] = gc_am_bm_estimate_1(x, prob_dist, param_pre, Y_k, mu, L, thr_val, count_break)
% GC_AM_BM_ESTIMATE used to estimate am and bm parameter, the result of
% this code obtains from second order equation of bm 
%
% INPUTS
% x             : the interval of x that is used to estimate parameters
% prob_dist     : distribution of different steps at given x
% param_pre.am  : am parameter of previous step
% param_pre.bm  : bm parameter of previous step
% Y_k           : observation
% mu            : mu parameter
% L             : eigenvector parameter
% thr_val     : threshold value for stop updaiting
% count_break : maximum number of updating 
%
% OUTPUTS
% Param_ab.am_estimate : estimated am for the next step
% Param_ab.bm_estimate : estimated bm for the next step



if nargin < 6
    error('the first 6 inputs "x, prob_dist, param_pre, Y_k, mu, L" should be passed')
end
if nargin == 6
    thr_val = 10^-3;
    count_break = 5;
end
if nargin == 7
    count_break = 5;
end

am_pre = param_pre.a;
bm_pre = param_pre.b;

update_val_b = zeros(length(bm_pre));
thr_exit = 1;
count = 0;
while thr_exit > thr_val && count < count_break
    count = count +1;
    
    % estimate am and bm by second order equation
    
    am_curr = am_estimator(bm_pre, am_pre, Y_k, mu, L, update_val_b, x, prob_dist);
%     am_curr = am_pre;
    
    [bm_curr, update_val_b] = bm_estimator(bm_pre, am_curr, Y_k, mu, L, x, prob_dist);
    
    
    % calculate difference between the current value and previous value
%     diff_am(count) = sum(abs(am_curr - am_pre));
    diff_bm(count) = sum(abs(bm_curr - bm_pre));
    
    % if the current and previous values are so close, change thr_exit in order to break
    if  diff_bm(count) < thr_val
        thr_exit = thr_val;
    end
    
    % update value of am and bm for the next loop
%     am_pre = am_curr;
    bm_pre = bm_curr;

end

% final results for am and bm
param_ab_estimate.a = am_curr;
param_ab_estimate.b = bm_curr;

am_curr;
bm_curr;

diff_am = param_pre.a - am_curr;
diff_bm = param_pre.b - bm_curr;

end

%%%%%%%%%%%%%%%%%%%%%%% nested fnctions %%%%%%%%%%%%%%%%%%%%%%%%%%

function [am_curr] = am_estimator(bm_pre, am_pre, Y_k, mu, L, update_val_b, x, prob_dist)
am_curr = am_pre;
K = length(Y_k);


% expection of x by given probablility distribution for all steps (K)
E_x = gc_mean_x(x, prob_dist); 

K = length(E_x); % number of all steps
T = sum(E_x(:)); % sum of all Exptected value of x for all steps

% estimation of C0_k, C1_k and C2_k (for caclulating Taylor series of Expection x of bm*x)
[C0, C1, C2] = gc_C0_C1_C2(x, bm_pre, prob_dist);


for m=1 : length(bm_pre)-1 % differnet channels
    
    temp_am = 0;
    for k=1: K % different steps
        z_k = Y_k(k, :).'- mu;
        c_k = z_k'*L;
        c_mk_2 = c_k(m)*c_k(m)';
        
        
        A0(k) = c_mk_2*C0(m,k);
        A1(k) = c_mk_2*C1(m,k);
        A2(k) = c_mk_2*C2(m,k);
        
    end
    
    am_curr(m,1) = log((sum(A0)-sum(A1)*update_val_b(m,1) + .5*sum(A2)*update_val_b(m,1)^2)./K);
end

end

function [bm_curr, update_val_b] = bm_estimator(bm_pre, am_pre, Y_k, mu, L, x, prob_dist)

% estimation of C0_k, C1_k and C2_k (for caclulating Taylor series of Expection x of bm*x)
[C0, C1, C2] = gc_C0_C1_C2(x, bm_pre, prob_dist);

E_x = gc_mean_x(x, prob_dist);
K = length(E_x); % number of all steps
T = sum(E_x(:)); % sum of all Exptected value of x for all steps

update_val_b = zeros(length(bm_pre),1);


for m=1 : length(bm_pre)-1 % differnet channels
    
    
    term_b1 = zeros(K,1);
    term_b2 = zeros(K,1);
    for k=1: K % different steps
        
        z_k = Y_k(k, :).'- mu;
        c_k = z_k'*L;
        c_mk_2 = c_k(m)*c_k(m)';
        
        A1(k) = c_mk_2*C1(m,k);
        A2(k) = c_mk_2*C2(m,k);
        
    end
    update_val_b(m,1) = (sum(A1) - exp(am_pre(m,1))*T)./sum(A2);
end
bm_curr = bm_pre + update_val_b;

end

%% C0, C1 and C2 calculator
function [C0, C1, C2] = gc_C0_C1_C2(x, bm, prob_dist)
% function gc_C0_C1_C2 calculate C0, C1 and C2 that is used to calculate am
% and bm
% 
% INPUTS     
% x          : interval of x that we want to evaluate it
% bm         : bm value 
% prob_dist  : estimated probability of x for all steps
%  
% OUTPUTS
% C0         : Expectaion of x of (exp(-bm*x))
% C1         : Expectaion of x of (x*exp(-bm*x))
% C2         : Expectaion of x of (x^2*exp(-bm*x))

if nargin < 3
    error('all inputs "x, bm, prob_dist" should be passed')
end

m = size(prob_dist);
% m(1) = length x
% m(2) = K (number of states)


C0 = zeros(length(bm) , m(2));
C1 = zeros(length(bm) , m(2));
C2= zeros(length(bm) , m(2));

for i=1:length(bm)
    f_x_C0 = exp(-bm(i).*x );
    f_x_C1 = x.*exp(-bm(i).*x );
    f_x_C2 = x.*x.*exp(-bm(i).*x );
    
    % Expected value of each element of D^-1
    E_lambda_C0_k = zeros(1, m(2));
    E_lambda_C1_k = zeros(1, m(2));
    E_lambda_C2_k = zeros(1, m(2));
    for j=1: m(2)
        E_lambda_C0_k(j) = f_x_C0'* prob_dist(:,j);
        E_lambda_C1_k(j) = f_x_C1'* prob_dist(:,j);
        E_lambda_C2_k(j) = f_x_C2'* prob_dist(:,j);
    end
    
    C0(i,:) = E_lambda_C0_k;
    C1(i,:) = E_lambda_C1_k;
    C2(i,:) = E_lambda_C2_k;
end


end

%% calculate expected value of x
function [expected_x] = gc_mean_x(x, prob_dist)
% GC_MEAN_X used to calculate Expected Value of x (mean)
%
% INPUTS
% x             : interval of x that we want to calculate expected value for it
% prob_dist     : given probablility distribution
%
% OUTPUTS
% expected_x   : expected value of x (mean) for given probability distribution

if nargin < 2
    error('all inputs "x, prob_dist" should be passed')
end

m = size(prob_dist);
% m(1) = length x
% m(2) = K (number of states)

expected_x = zeros(m(2),1); % mean
for i=1: m(2)
    expected_x(i) = x'*prob_dist(:,i);
end

end

