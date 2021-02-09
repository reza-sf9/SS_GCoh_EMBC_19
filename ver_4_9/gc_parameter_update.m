function [Param_update] = gc_parameter_update(Param, ParamaUpdate, Bayes)
% GC_PARAMETER_UPDATE this function is used to estimate parameters
% 
% INPUTS
% Param                        : is a structure that contains parameters
%    Param.sv                  : variance 
%    Param.s_ab                : a and b of random walk model
%    Param.mu                  : mean
%    Param.L                   : eigenvalue
%    Param.ab                  : a structure of am and bm
%          Param.ab.a          : am
%          Param.ab.b          : bm
% ParamaUpdate                 : is a structure specifies which parameters should be updated
%     ParamaUpdate.sv_update   : if is 1 sv updates
%     ParamaUpdate.mu_update   : if is 1 mu updates
%     ParamaUpdate.L_update    : if is 1 L updates
%     ParamaUpdate.ab_update   : if is 1 am and bm updates
% Bayes                        : a structure of bayes filter-smoother results
%     Bayes.smoother           : smoother result 
%     Bayes.filter             : filter result 
%     Bayes.oneStep            : one step prediction of fiter
%     Bayes.filter_prior       : prior probability of filter 
%     Bayes.stateTransient     : state transient (step k to k+1)t 
%     Bayes.x                  : an interval of x that we want to evalueate
%     Bayes.Y_k                : observation (measured of FFT)
% 
% OUTPUTS
% Param_update                        : is a structure contains updated values
%         Param_update.s_ab           : updated value of a and b parameters (of random walk model)
%         Param_update.sv             : updated value of variance
%         Param_update.mu             : updated value of mu
%         Param_update.L              : updated value of L (eigenvlaues)
%         Param_update.ab             : an sutructure contains updated value of am and bm 
%                 Param_update.ab.a   : updated value of am
%                 Param_update.ab.b   : updated value of bm

if nargin < 3
    error('the first 3 inputs "ParamaUpdate, Param, Bayes" should be passed')
end

% update of don't update info
sv_update = ParamaUpdate.sv_update;

sv = Param.sv;

% structure of bayes filter-smoother results
smoother_estimate = Bayes.smoother;    % result of smoother
filter_estimate = Bayes.filter;        % result of filter
f_xx = Bayes.stateTransient;           % state transient (step k to k+1)
oneStep_prediction = Bayes.oneStep;    % oneStep prediction of filter
x = Bayes.x;                           % an interval of x that we want to evalueate
p0 = Bayes.filter_prior;               % prior probability is used in filter estimation
    
if sv_update == 1
    % estimate sv2
    sv = gc_sv2_estimate(x, f_xx, filter_estimate,  smoother_estimate, oneStep_prediction, p0);
end

% update mu, L and am & bm 
thr_val = 1e-6;
count_break = 10;
[mu_estimate, L_estimate, param_ab_estimate] = gc_update_mu_L_ab(ParamaUpdate, Param, Bayes, thr_val, count_break);

% update 
Param_update.s_ab = Param.s_ab;
Param_update.sv = sv;
Param_update.mu = mu_estimate;
Param_update.L = L_estimate;
Param_update.ab = param_ab_estimate;

end