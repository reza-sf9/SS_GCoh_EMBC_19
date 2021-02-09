function [mu, L, param_ab] = gc_update_mu_L_ab(ParamaUpdate, Param, Bayes, thr_val, count_break)
% GC_UPDATE_MU_L_AB estimates mu, L, am and bm parameters
% 
% INPUTS
% ParamaUpdate                 : is a structure specifies which parameters should be updated
%     ParamaUpdate.sv_update   : if is 1 sv updates
%     ParamaUpdate.mu_update   : if is 1 mu updates
%     ParamaUpdate.L_update    : if is 1 L updates
%     ParamaUpdate.ab_update   : if is 1 am and bm updates
% Param                        : is a structure that contains parameters
%    Param.sv                  : variance 
%    Param.s_ab                : a and b of random walk model
%    Param.mu                  : mean
%    Param.L                   : eigenvalue
%    Param.ab                  : a structure of am and bm
%          Param.ab.a          : am
%          Param.ab.b          : bm
% Bayes                        : a structure of bayes filter-smoother results
%     Bayes.smoother           : smoother result 
%     Bayes.filter             : filter result 
%     Bayes.oneStep            : one step prediction of fiter
%     Bayes.filter_prior       : prior probability of filter 
%     Bayes.stateTransient     : state transient (step k to k+1)t 
%     Bayes.x                  : an interval of x that we want to evalueate
%     Bayes.Y_k                : observation (measured of FFT)
% thr_val                      : a threshold value for stopping loop
% count_break                  : maximum nmber of updating loop
% 
% OUTPUTS
% mu                           : estimated mu
% L                            : estimated L (eigenvector)
%    Param.ab                  : a structure of estimated values of  am and bm
%          Param.ab.a          : estimated value of am
%          Param.ab.b          : estimated value of bm

if nargin < 3
    error('the first 3 inputs "ParamaUpdate, Param, Bayes" should be passed')
end
if nargin == 3
    thr_val = 0.01;
    count_break = 5;
end
if nargin == 4
    thr_val = 0.01;
    count_break = 5;
end

% extracting from structures
mu_update = ParamaUpdate.mu_update;
L_update = ParamaUpdate.L_update;
ab_update = ParamaUpdate.ab_update;

mu = Param.mu;
L = Param.L ;
param_ab.a = Param.ab(:,1);
param_ab.b = Param.ab(:,2);

smoother_estimate = Bayes.smoother; % result of smoother
x = Bayes.x;                        % an interval of x that we want to evalueate
Y_k = Bayes.Y_k;                    % observation (FFT measure)

error_mu = 0;
error_L = 0;
error_ab = 0;

count = 0;
thr = 100;
while thr > thr_val && count < count_break
    count = count+1;
    
    if mu_update == 1
        % estimate mu
        [mu, mu_diff] = gc_mu_estimate(x, param_ab, smoother_estimate, L, Y_k, mu);
        error_mu(count) = abs(sum(mu_diff(:)));
    end
    
    
    if L_update == 1
        % estimate L
        [L, L_diff] = gc_L_estimate(x, param_ab, smoother_estimate, L, Y_k, mu,1e-6,10);
        error_L(count) = abs(sum(sum(L_diff)));
    end
    
    if ab_update==1
        % estimate am & bm
        
        % second order solve
%         [param_ab, am_diff, bm_diff] = gc_am_bm_estimate(x, smoother_estimate, param_ab, Y_k, mu, L);
        
        % first order solve
        [param_ab, am_diff, bm_diff] = gc_am_bm_estimate_1(x, smoother_estimate, param_ab, Y_k, mu, L);
        error_ab(count) = abs(sum(am_diff)) + abs(sum(bm_diff));
    end
    
    error_tot(count) = error_mu(end) + error_L(end) + error_ab(end);
    thr = error_tot(count);
    
end
param_ab = [param_ab.a param_ab.b];

end