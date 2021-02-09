function [smoother_estimate]= gc_bayes_smoother(f_xx, oneStep_prediction, filter_estimate)
% GC_BAYES_SMOOTHER used to estimate the posterior probability based Bayes
% Smoother
% (for further information please take a look at part 1-6 of supplemntary)
%
% INPUTS
% f_xx_smoother     : state transient matrix
% oneStep_prediction    : one step prediction of Bayes Filter
% filter_estimate  : posterior probability of Bayes Filter
%
% OUTPUTS
% posterior_smoother: posterior probability of Bayes Smoother

 
if nargin < 3
    error('the first 3 inputs "f_xx_smoother, oneStep_prediction, filter_estimate" should be passed')
end


m = size(filter_estimate);
% m(1) = length(x)
% m(2) = K (number of states)
K = m(2);

L_K_smoother = zeros(m(1), m(2));     % Likelihood
oneStep_smoother = zeros(m(1), m(2)); % One-Step Prediction
smoother_estimate = zeros(m(1), m(2)); % One-Step Prediction

smoother_estimate(:,end) = filter_estimate(:,end);

% solve problem of having NAN results
epsilon_v = 10^-9;
oneStep_prediction = max(oneStep_prediction ,epsilon_v); 

for i=K-1:-1: 1
    
    % posterior of early step
    posterior_early = smoother_estimate(:,i+1);
    
    % Posterior
    posterior_current = (f_xx'*(posterior_early./oneStep_prediction(:,i+1))).* filter_estimate(:, i);  % eq (13) and eq(14)
    
    % normalizing
    posterior_current = posterior_current./ sum(posterior_current);
    
    smoother_estimate(:,i) = posterior_current;
    
end


end