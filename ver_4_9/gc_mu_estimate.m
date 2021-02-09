function [mu_estimate, diff_mu] = gc_mu_estimate(x, param, prob_dist, L, Y_k, mu_init)
% GC_MU_ESTIMATE estimates the mu value 
%
% INPUTS
% x         : an interval of all possible state variables 
% param     : used param (am, bm) for generatig lambda
% prob_dist : given probablility distribution
% L         : eigenvector matrix
% Y_k       : observation 
%
% OUTPUTS
% mu_estimate : estimated mu

if nargin < 5
    error('all inputs "x, param, prob_dist, L, Y_k" should be passed')
end

% get expected of D^-1
[E_D_inv] = gc_expectation_D_inv(x, param, prob_dist);


m = size(Y_k);
% m(1) = number of steps (K)
% m(2) = number of channels 

sum_1 = zeros(m(2), m(2));
sum_2 = zeros(m(2), 1);
for k=1 : m(1)
    sum_1 = sum_1 + L*E_D_inv(:,:,k)*L';
    
    sum_2 = sum_2 + (L*E_D_inv(:,:,k)*L')*Y_k(k,:).';
end

mu_estimate = inv(sum_1)*sum_2;

diff_mu = mu_estimate - mu_init;

end

%%%%%%%%%%%%%%%%%%%%%%%%%% Nested functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate expectation of D^-1
function [E_D_inv] = gc_expectation_D_inv(x, param, prob_dist)
% GC_EXPECTED_D_INV calculates expected value of D^-1, is used to estimate
% mu
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