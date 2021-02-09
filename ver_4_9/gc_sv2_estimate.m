function [sigma_2] = gc_sv2_estimate(x, f_xx, filter_estimate,  smoother_estimate, oneStep_prediction, p0)
% GC_SV2_ESTIMATE is used to estimate variance
%
% INPUT
% var_x  : expected value of x^2
% corr_x : correlation between x(k) and x(k+1)
%
% OUTPUT
% sigma_2 = estimated variance

if nargin < 3
    error('all inputs "var_x,var_0, corr_xt" should be passed')
end

var_0 = (x.*x)'*p0;
var_x = gc_var_x(x, smoother_estimate);
joint_prob = gc_joint_prob(f_xx, filter_estimate,  smoother_estimate, oneStep_prediction, p0);
corr_x = gc_cov_x(x, joint_prob);

K = length(var_x);

temp = 2*sum(var_x)+var_0-var_x(end)-2*sum(corr_x);
sigma_2 = temp/K;
sigma_2 = max(1e-6,sigma_2);


end

%%%%%%%%%%%%%%%%%%%%%%%%%% Nested functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate var of x
function [expected_x2] = gc_var_x(x, prob_dist)
% GC_var_X used to calculate Expected Value of x^2 (var)
% 
% INPUTS
% x                   : interval of x that we want to calculate expected value for it
% prob_dist           : given probablility distribution
% 
% OUTPUTS
% expected_x2   : expected value of x^2 (var) for given probability distribution

if nargin < 2
    error('all inputs "x, prob_dis" should be passed')
end

m = size(prob_dist);
% m(1) = length x
% m(2) = K (number of states)
xx = x.*x;
expected_x2 = zeros(m(2),1);% var
for i=1: m(2)
    expected_x2(i) = xx'* prob_dist(:,i);
end


end

%% calculate joint probability distribution of step k and k+1
function [joint_prob] = gc_joint_prob(f_xx, filter_estimate,  smoother_estimate, oneStep_prediction, p0)
% GC_JOINT_PROB used for calculate joint probability of x(k) and x(k+1)
% 
% INPUTS
% transient           : transient matrix p(x(k+1)|x(k))
% filter_estimate     : posterior probability of filter
% smoother_estimate   : posterior probability of smoother
% oneStep_prediction  : one Step prediction of filter
% p0                  : prior probablitly
%
% OUTPUTS
% p_joint  : joint porobability of step k and k+1

if nargin < 5
    error('all inputs "f_xx, filter_estimate,  smoother_estimate, oneStep_prediction, p0" should be passed')
end

epsilon_v = 10^-9;
oneStep_prediction = max(oneStep_prediction ,epsilon_v);


m = size(smoother_estimate);


joint_prob = zeros(m(1), m(1), m(2));

joint_prob(:,:,1) = repmat((smoother_estimate(:,1)./ oneStep_prediction(:,1)),1,m(1)) .* f_xx .* repmat(p0',m(1),1);
joint_prob(:,:,1) = joint_prob(:,:,1)/sum(sum(joint_prob(:,:,1)));

for k=1: m(2)-1
    joint_prob(:,:,k+1)=  repmat((smoother_estimate(:,k+1)./ oneStep_prediction(:,k+1)),1,m(1)) .* f_xx .* repmat(filter_estimate(:,k)',m(1),1);
    joint_prob(:,:,k+1) = joint_prob(:,:,k+1)/sum(sum(joint_prob(:,:,k+1)));
end


end

%% caclulate correlation between step k and k+1
function cov_x = gc_cov_x(x,joint_dist)
% GC_cov_X used to calculate Covariance matrix with using joint porbability
% 
% INPUTS
% x           : interval of x that we want to calculate expected value for it
% joint_dist  : given joint probablility distribution
% 
% OUTPUTS
% cov_x       : covariance matrix of x for given joint probability distribution

if nargin < 2
    error('all inputs "x,joint_dist" should be passed')
end

m = size(joint_dist);
% m(1) = length x
% m(2) = length x
% m(3) = K (number of states)
cov_x  = zeros(m(3),1);

for i=1:m(3)
    joint_prob_temp = joint_dist(:,:,i);
    cov_x(i) = sum(sum(repmat(x,1,m(1)).*joint_prob_temp.*repmat(x',m(1),1)));
end

end
