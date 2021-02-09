function [filter_estimate, oneStep_prediction ,L_K_filter, prior] = gc_bayes_filter(f_xx, Y_k, x, L, mu, param)
% GC_BAYES_FILTER this function estimate the posterior probability of the
% state by estimating the Likelihood and One-Step Prediction distributions
% (for further information please take a look at part 1-5 of supplemntary)
%
% INPUTS
% f_xx      : state transient matrix
% p_xx_pre  : distribution of previous posterior
% y_k       : the current observation
% x         : an interval of state variables
% L         : given eigenvector
% mu        : given mu
% param.a   : a paramter for generating lambda
% param.b   : b paramter for generating lambda
%
% OUTPUTS
% filter_estimate    : distribution of current posterior
% oneStep_prediction : one step prediction
% L_K_filter         : likelihood
% prior              : prior priorirty used for estimating filter

if nargin < 6
    error('all inputs "f_xx, Y_k, x, L, mu, param" should be passed')
end

%%

m = size(Y_k);
% m(1) = K (number of states)
% m(2) = number of channels
K = m(1);

prior = gc_prior_filter(x, 2);


L_K_filter      = zeros(length(x), K);     % Likelihood
oneStep_prediction  = zeros(length(x), K); % One-Step Prediction

filter_estimate  = zeros(length(x),K);
% estimation of posterior probablity of x based BAYES FILTER
for i=1:K
    y_k = Y_k(i,:).';  % current observation
    
    if i==1
       % call One-Step Prediction calculator
       oneStep = gc_step_predict(f_xx, prior);
    else
       % call One-Step Prediction calculator
       oneStep = gc_step_predict(f_xx, filter_estimate(:,i-1)); 
    end
    
    % call Likelihood calculator
    L_k = gc_likelihood(y_k, x, L, mu, param);
    
    % Posterior distribution
    posterior_curr = oneStep.*L_k;           % eq (9)
    
    % normalizing Posterior
    posterior_curr = posterior_curr./sum(posterior_curr.');
    
    % normalized L_k
    L_k = L_k/max(eps,sum(L_k));
     

    filter_estimate(:,i) = posterior_curr; % Posterior
    L_K_filter(:,i)       = L_k;               % Likelihood
    oneStep_prediction(:,i)   = oneStep;           % One-Step Prediction
    %% in index i, i keep p(x_i|1:i-1)
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%% Nested functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% one step prediction
function [f_xx_1] = gc_step_predict(f_xx, p_x)
% GC_STEP_PREDICT calculates posterior probability
% (for further information please take a look at part 1-5-1 of supplemntary)
% 
% INPUTS
% f_xx     : state transient matrix
% p_x      : posterior of the previous step
%
% OUTPUTS
% f_xx_1  : one step prediction of Bayes filter

if nargin < 2
    error('all inputs "f_xx and p_x" should be passed')
end

% calculate one step prediction from the last posterior and state transient
f_xx_1 = f_xx*p_x;           % eq (10)

end

%% Likelihood claculator
function [L_k] = gc_likelihood(y_k, x, L, mu, param)
% GC_LIKELIHODD calculates likelihood distribution for each given
% observatoin (y_k) in an interval of state variables (x)
% (for further information please take a look at part 1-5-2 of supplemntary)
% 
% INPUTS
% y_k      : the current observation
% x        : an interval of state variables
% L        : given eigenvector
% mu       : given mu
% param.a  : a paramter for generating lambda 
% param.b  : b paramter for generating lambda 
% 
% OUTPUTS
% L_k    : Likelihood distribution(a vector) 

if nargin < 5
    error('all inputs "y_k, x, L, mu and param" should be passed')
end



% Likelihood 
L_k = zeros(length(x), 1);

for i=1:length(x) % for differnt states
    
    % eigenvalues matrix for x_k
    D_xk = diag(exp(param.a + x(i).*param.b));
    
    % covriance matrix for x_k
    Gamma = L*D_xk*L';
    
    % calc likelihood for given observation and state
    L_k(i,1) = exp(-real((y_k-mu)'*pinv(Gamma)*(y_k-mu)))./det(D_xk);   % eq (11)
end



end

%% prior probability generator
function prior = gc_prior_filter(x, opption_prior,mu, sigma)
% GC_PRIOR_FILTER generate prior probability of filter
%
% INPUTS
% x            : an interval of state variables
% option_prior : % we have 2 options  for prior probability, first is a uniform prior prob and
%                 second is based a normal distribution 
% mu           : given mu
% sigma        : given variance
%
% OUTPUTS
% prior        : generated probability for prior step of filter


if nargin == 0
    error('the first  inputs "x" should be passed')
end
if nargin == 1
    opption_prior =1;
    mu = 0;
    sigma = sqrt(0.1);
end
if nargin == 2
    mu = 0;
    sigma = sqrt(0.1);
end
if nargin == 3
    sigma = sqrt(0.1);
end

switch opption_prior
    case 1
        % 1- uniform prior
        prior = ones(length(x),1)/length(x);
    case 2
        % 2- using a normal dist for prior probability
        prior = normpdf(x, mu, sigma);
        prior = prior./sum(prior);
end

end
