function [prior, transmat, mixmat, mu, sigma] = hmm_em(sm, tracks, varargin)
%
%does hmm parameter estimation for the segmentation model
%function [prior, transmat, mixmat, mu, sigma] = hmm_em(sm, tracks,
%varargin)
%
%optional args:
%   'update', true/false (default true); if true, assigns new values to sm



update = true;
adj_prior = false;
varargin = assignApplicable(varargin);

%{
if isempty(sm.hmm_musub) || isempty(sm.hmm_sigmnorm)
    warning('you will probably have better results if you standardize data set');
    stand = false;
    mu = [];
    sigma2 = [];
else
    stand = true;
    mu = sm.hmm_musub;
    sigma2 = sm.hmm_sigmnorm;
end

[data, known_points] = sm.createHMMData(tracks, 'stand', stand, 'mu', mu, 'sigma2', sigma2);
%}
[data, known_points] = sm.createHMMData(tracks); %standardization is taken care of with createHMMData
[ll_trace, prior, transmat, mu, sigma, mixmat] = mhmm_em_gershow(data, ...
    sm.hmm_prior, sm.hmm_transmat, sm.hmm_mean, sm.hmm_cov, sm.hmm_mix, ...
    'adj_trans', logical(sm.allowedTransitions), 'known_points', known_points, 'adj_prior', adj_prior, varargin{:});

%{
[ll_trace, prior, transmat, mu, sigma, mixmat] = mhmm_em(data, ...
    sm.hmm_prior, sm.hmm_transmat, sm.hmm_mean, sm.hmm_cov, sm.hmm_mix, 'max_iter', 10, 'adj_prior', false) ;
%}
plot (ll_trace);

if (update)
    sm.hmm_prior = prior;
    sm.hmm_transmat = transmat;
    sm.hmm_mean = mu;
    sm.hmm_cov = sigma;
    sm.hmm_mix = mixmat;
end