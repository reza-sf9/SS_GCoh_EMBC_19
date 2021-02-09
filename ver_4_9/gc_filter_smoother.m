function [Bayes, pmf_gc, s] = gc_filter_smoother(Y_k, x, Param)
% GC_FILTER_SMOOTHER estimates posterior with bayes filter and smoother ,
% after calculating smoother, we estimate a te distribution for GC
%
% INPUTS
% Y_k                 : observation (measured of FFT)
% x_k                 : state variables
% Param               : is a structure that contains parameters
%    Param.sv         : variance 
%    Param.s_ab       : a and b of random walk model
%    Param.mu         : mean
%    Param.L          : eigenvalue
%    Param.ab         : a structure of am and bm
%          Param.ab.a : am
%          Param.ab.b : bm
% 
% OUTPUTS
% Bayes                    : a structure of bayes filter-smoother results
%     Bayes.smoother       : smoother result 
%     Bayes.filter         : filter result 
%     Bayes.oneStep        : one step prediction of fiter
%     Bayes.filter_prior   : prior probability of filter 
%     Bayes.stateTransient : state transient (step k to k+1)t 
%     Bayes.x              : an interval of x that we want to evalueate
%     Bayes.Y_k            : observation (measured of FFT)
% pmf_gc                   : distribution of Global Coherence
% s                        : an interval for calculating GC


% state transient probabilities - p(x(k)|x(k-1))
a = Param.s_ab(1);
b = Param.s_ab(2);
sv = Param.sv;
mu = Param.mu;
L = Param.L ;
param_ab.a = Param.ab(:,1);
param_ab.b = Param.ab(:,2);

% state transient probabilities - p(x(k)|x(k-1))
f_xx = gc_state_transition(x, sv ,a ,b);

% Bayes filter
[filter_estimate, oneStep_prediction, L_K_filter, p0] = gc_bayes_filter(f_xx, Y_k, x, L, mu, param_ab);

% Bayes smoother
[smoother_estimate]= gc_bayes_smoother(f_xx, (oneStep_prediction), (filter_estimate));

% pmf distribution of Global Coherence
[pmf_gc, s] = gc_pdf_gc(x, mu, param_ab, smoother_estimate);

% structure of bayes filter-smoother results
Bayes.smoother = (smoother_estimate);
Bayes.filter = filter_estimate;
Bayes.oneStep = oneStep_prediction;
Bayes.filter_prior = p0;
Bayes.stateTransient = f_xx;
Bayes.x = x;
Bayes.Y_k = Y_k;


end