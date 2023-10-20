function initHMM (sm, varargin)
%make a bayesnet hidden markov model based on segmentation model
%function bnet = createBNET (sm, varargin)

%this is the topography of a mixture of gaussians hidden markov model
%node 1 = Q at time t informs node 1 = Q at time t+1;  this is the
%only inter-slice connection
%
%node 1 at Q informs node 2 - the mixture ratio, and node 3, the observed
%node 2 informs node 3;  these are the intra-slice connections

numgaussians = 1;
nglist = [];
tracklist = [];
varargin = assignApplicable(varargin);

if (isempty(sm.hmm_musub) || isempty(sm.hmm_sigmnorm))
    setStandardization(sm, tracklist);
end

Q = length(sm.segmentationClusters); % num hidden states
O = length(sm.segmentationClusters(1).datafields); % size of observed vector
M = numgaussians; % num mixture components per state

prior = [sm.segmentationClusters.priorProbability];
if (length(prior) ~= length(sm.segmentationClusters) || any(~isfinite(prior)))
    prior = zeros(size(sm.segmentationClusters));
    for j = 1:length(sm.segmentationClusters)
        prior(j) = length([sm.segmentationClusters(j).knownPoints.inds]);
    end
end
%prior
prior = prior./sum(prior); %make sure prior probabilities sum to 1
for j = 1:length(prior)
    sm.segmentationClusters(j).priorProbability = prior(j);
end
sm.hmm_prior = prior;

args = {'stand', true, 'mu', sm.hmm_musub, 'sigma2', sm.hmm_sigmnorm};

mixmat = zeros([Q M]);
mu = zeros(O, Q, M);
sigma = zeros(O, O, Q, M);
[cm, ccv, cmr] = meanAndCovOfKnownPoints(sm.segmentationClusters, M, args{:});

for j = 1:length(sm.segmentationClusters)
    mixmat(j,:) = cmr{j};
    mu(:,j,:) = cm{j}';
    sigma(:,:,j,:) = ccv{j};
end
sm.hmm_mix = mixmat;
sm.hmm_mean = mu;
sm.hmm_cov = sigma;

transmat = zeros(Q,Q);

for i = 1:Q
    transmat(i,:) = prior .* sm.allowedTransitions(i,:) ./ (sum(prior.*sm.allowedTransitions(i,:)));
end
sm.hmm_transmat = transmat;  
%transmat
sc = sm.segmentationClusters;

%improve the transmat & prior guess based on the track with the most points defined
kp = [sc.knownPoints];
tr = unique([kp.track]);
nump = zeros(size(tr));
for j = 1:length(tr)
    nump(j) = length([kp([kp.track] == tr(j)).inds]);
end
[~,I] = max(nump);
tr = tr(I);
%I

%adjust the transmat
disp('fitting transmat');
sm.hmm_em(tr, 'adj_prior', 0, 'adj_mix', 0,'adj_mu', 0, 'adj_Sigma',0,'adj_trans',1);
sm.hmm_transmat
%adjust the prior
disp('fitting prior');
sm.hmm_em(tr, 'adj_prior', 1, 'adj_mix', 0,'adj_mu', 0, 'adj_Sigma',0,'adj_trans',0);
sm.hmm_transmat
if (isempty(tracklist))
    return;
end
if (isempty(nglist))
    nglist = repmat(M, size(sc));
end
%now go through all tracks in tracklist and assign each point to a best
%guess based on viterbi

%fit all points to a first guess 
gs = cell(size(tracklist));
for j = 1:length(tracklist)
    gs{j} = sm.doViterbi(tracklist(j));
    sum(gs{j} == 2)
end

%cluster the guessed points into nglist(j) gaussians
sc2 = SegmentationCluster();
sc2.datafields = sc(1).datafields;
sc2.operation = sc(1).operation;
clear kp;
for j = 1:length(gs)
    kp(j).track = tracklist(j);
end 
cm = cell(size(sc));
ccv = cell(size(sc));
cmr = cell(size(sc));
for k = 1:length(sc)
    sc(k).name
    for j = 1:length(gs)
        kp(j).inds = find(gs{j} == k);
    end
    sc2.knownPoints = kp;
    
    try
        [cm{k}, ccv{k}, cmr{k}] = meanAndCovOfKnownPoints(sc2, nglist(k), args{:});
    catch
        disp(['failed to cluster ' sc(k).name ' properly with ' num2str(nglist(k)) ' clusters -- using single cluster instead']);
        [cm{k}, ccv{k}, cmr{k}] = meanAndCovOfKnownPoints(sc2, 1, args{:});
    end
end    

%initialize all to zero, since we may have different number of clusters
%then copy in the fitted clusters to the overall model
M = max(nglist);
mixmat = zeros([Q M]);
mu = zeros(O, Q, M);
sigma = ones(O, O, Q, M);
for j = 1:length(sm.segmentationClusters)
    inds = 1:length(cmr{j});
    mixmat(j,inds) = cmr{j};
    %size(mu(:,j,inds))
    %size(cm{j})
    mu(:,j,inds) = cm{j}';
    sigma(:,:,j,inds) = ccv{j};
end
sm.hmm_mix = mixmat;
sm.hmm_mean = mu;
sm.hmm_cov = sigma;
