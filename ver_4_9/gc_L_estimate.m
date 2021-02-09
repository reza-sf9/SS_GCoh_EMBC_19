function [L_estimate, L_diff] = gc_L_estimate(x, param, prob_dist,L ,Y_k, mu, thr_val, count_break)
% GC_L_ESTIMATE estimates the eigenvector matrix (L)
%
% INPUTS
% x           : an interval of all possible state variables
% param       : used param (am, bm) for generatig lambda
% prob_dist   : given probablility distribution
% L           : eigenvector matrix
% Y_k         : observation
% mu          : mean of model
% thr_val     : threshold value for stop updaiting
% count_break : maximum number of updating
%
% OUTPUTS
% L_estimate : estimated mu

if nargin < 6
    error('all inputs "x, param, prob_dist,L ,Y_k, mu" should be passed')
end
if nargin == 6
    thr_val = 10^-2;
    count_break = 5;
end
if nargin == 7
    count_break = 5;
end

L_init = L;

% get expected of D^-1
[E_D_inv] = gc_expectation_D_inv(x, param, prob_dist);

m = size(E_D_inv);
% m(1), m(2) = number of channels
% m(3) = number of state variables (K)

I = eye(m(2));

L_pre = L;
thr_L = 100;
count = 0;
beta_1 = 50;
beta_2 = 50;
eps_u = -10^-5;
eps_v = eps_u;

while (thr_L > thr_val) && (count < count_break)
    count = count+1;
    
    U_pre = real(L_pre);
    V_pre = imag(L_pre);
    
    p=1;
    [norm_2_p_gamma, norm_2_p_theta] = gc_getting_norm(U_pre,V_pre,I,p);
    
    % update U
    U_curr = gc_update_U(U_pre, V_pre, Y_k, mu, E_D_inv, norm_2_p_gamma, norm_2_p_theta, beta_1, beta_2, eps_u);
    
    % update V
    V_curr = gc_update_V(U_curr, V_pre, Y_k, mu, E_D_inv, norm_2_p_gamma, norm_2_p_theta , beta_1, beta_2, eps_v);
    
    L_curr = U_curr + V_curr*1i;
    
    Gamma = U_curr*U_curr'+V_curr*V_curr'-I
    Thet = U_curr*V_curr'-V_curr*U_curr'

    L_curr*L_curr'
    
    diff_L = abs(L_curr-L_pre);
    thr_L = real(sum(diff_L(:))) + imag(sum(diff_L(:)));
    L_pre = L_curr;
end

L_estimate = L_curr;
L_diff = L_estimate - L_init;
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Nested functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculated expectation of D^-1
function [E_D_inv] = gc_expectation_D_inv(x, param, prob_dist)
% GC_EXPECTED_D_INV calculates expected value of D^-1, is used to estimate mu
%
% INPUTS
% x         : an interval of all possible state variables
% param     : used param (am, bm) for generatig lambda
% prob_dist : given probablility distribution
%
% OUTPUTS
% E_D_inv   : expected value of D^-1 for all steps

if nargin < 3
    error('all inputs "x, param, prob_dist" should be passed')
end

m = size(prob_dist);
% m(1) = length x
% m(2) = K (number of states)


E_lambda_inv = zeros(length(param.a) , m(2));
% row = differnt am and bm (length = num of channels)
% col = differnt steps (length = K)

for i=1:length(param.a)
    f_x = exp(-param.a(i) -param.b(i).*x );
    
    % Expected value of each element of D^-1
    E_lambda_inv_temp = zeros(1, m(2));
    for j=1: m(2)
        E_lambda_inv_temp(j) = f_x'* prob_dist(:,j);
    end
    
    E_lambda_inv(i,:) = E_lambda_inv_temp;
end


E_D_inv = zeros(length(param.a),length(param.a),m(2));
for i=1:m(2)
    E_D_inv(:,:,i) = diag(E_lambda_inv(:,i));
end


end

function [U_curr] = gc_update_U(U, V, Y_k, mu, E_D_inv, norm_Gamma, norm_Theta, beta_1, beta_2, eps_u)

if nargin < 5
    error('all inputs "U, V, Y_k, mu, E_D_inv" should be passed')
end
if nargin == 5
    beta_1 = 1;
    beta_2 = 1;
    eps_u = 10^-2;
end
if nargin == 6
    beta_2 = 1;
    eps_u = 10^-2;
end
if nargin == 7
    eps_u = 10^-2;
end

m = size(E_D_inv);

I = eye(m(1));
dev_R = zeros(m(1), m(1));

for k = 1: m(3)
    
    r_k = real(Y_k(k, :).'- mu);
    i_k = imag(Y_k(k, :).'- mu);
    C_k = E_D_inv(:,:,k);
    
    dev_R = dev_R + (r_k*r_k'+ i_k*i_k')*U*C_k + (r_k*i_k'- i_k*r_k')*V*C_k;
    
end

dev_Gamma = 4*(U*U'+V*V'-I)*U;
dev_Theta = 4*(U*V'-V*U')*V;

dev_Gamma_p = dev_Gamma./norm_Gamma;
dev_Theta_p = dev_Theta./norm_Theta;

% beta_1 = max(2*dev_R);
% beta_2 = max(2*dev_R);
k = .1;
beta_1 = k*(max(2*dev_R)./max(dev_Gamma_p));
beta_2 = k*(max(2*dev_R)./max(dev_Theta_p));

% grad_val = 2*dev_R+ beta_1*dev_Gamma + beta_2*dev_Theta;
grad_val = 2*dev_R+ beta_1*dev_Gamma_p + beta_2*dev_Theta_p;


U_curr = U - eps_u*grad_val;

end

function [V_curr] = gc_update_V(U, V, Y_k, mu, E_D_inv, norm_Gamma, norm_Theta, beta_1, beta_2, eps_v)

if nargin < 5
    error('all inputs "U, V, Y_k, mu, E_D_inv" should be passed')
end
if nargin == 5
    beta_1 = 1;
    beta_2 = 1;
    eps_u = 10^-2;
end
if nargin == 6
    beta_2 = 1;
    eps_u = 10^-2;
end
if nargin == 7
    eps_v = 10^-2;
end

m = size(E_D_inv);

I = eye(m(1));
dev_R = zeros(m(1), m(1));
dev_Gamma = zeros(m(1), m(1));
dev_Theta = zeros(m(1), m(1));

for k = 1: m(3)
    
    r_k = real(Y_k(k, :).'- mu);
    i_k = imag(Y_k(k, :).'- mu);
    C_k = E_D_inv(:,:,k);
    
    dev_R = dev_R + (r_k*r_k'+ i_k*i_k')*V*C_k + (i_k*r_k'-r_k*i_k')*U*C_k;
    
end

dev_Gamma = 4*(V*V'+U*U'-I)*U;
dev_Theta = 4*(V*U'-U*V')*U;

dev_Gamma_p = dev_Gamma./norm_Gamma;
dev_Theta_p = dev_Theta./norm_Theta;

% beta_1 = max(2*dev_R);
% beta_2 = max(2*dev_R);
k = .1;
beta_1 = k*(max(2*dev_R)./max(dev_Gamma_p));
beta_2 = k*(max(2*dev_R)./max(dev_Theta_p));

% grad_val = 2*dev_R+ beta_1*dev_Gamma + beta_2*dev_Theta;
grad_val = 2*dev_R+ beta_1*dev_Gamma_p + beta_2*dev_Theta_p;

V_curr = V - eps_v*grad_val;

end


function [norm_2_p_gamma, norm_2_p_theta] = gc_getting_norm(U,V,I,p)

    % norm p 
    GAMMA = V*V'+U*U'-I;
    THETA = V*U'-U*V';
    
   
    norm_2_p_gamma = 0;% norm 2-p gamma
    norm_2_p_theta = 0;% norm 2-p theta
    
    q =2-p;
    for i=1:length(GAMMA)
        for j=1:length(THETA)
            norm_2_p_gamma = norm_2_p_gamma + abs(GAMMA(i,j))^q;
            norm_2_p_theta = norm_2_p_theta + abs(THETA(i,j))^q;
        end
    end

    norm_2_p_gamma = norm_2_p_gamma^(1/q);
    norm_2_p_theta = norm_2_p_theta^(1/q);
end


