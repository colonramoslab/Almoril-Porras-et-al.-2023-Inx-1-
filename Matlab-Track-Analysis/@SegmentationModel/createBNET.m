function bnet = createBNET (sm, varargin)
%make a bayesnet hidden markov model based on segmentation model
%function bnet = createBNET (sm, varargin)

%this is the topography of a mixture of gaussians hidden markov model
%node 1 = Q at time t informs node 1 = Q at time t+1;  this is the
%only inter-slice connection
%
%node 1 at Q informs node 2 - the mixture ratio, and node 3, the observed
%node 2 informs node 3;  these are the intra-slice connections

numgaussians = 1;
varargin = assignApplicable(varargin);

intra = zeros(3);
intra(1,[2 3]) = 1;
intra(2,3) = 1;
inter = zeros(3);
inter(1,1) = 1;



Q = length(sm.segmentationClusters); % num hidden states
O = length(sm.segmentationClusters(1).datafields); % size of observed vector
M = numgaussians; % num mixture components per state
Q
O
M
ns = [Q M O];
dnodes = [1 2]; %discrete nodes
onodes = 3; %observed nodes
eclass1 = [1 2 3]; %equivalence classes for first slice
eclass2 = [4 2 3]; %equivalence classes for all slices
%Q1 is different from Q2....N because Q1 does not have a prior state

bnet = mk_dbn(intra, inter, ns, 'discrete', dnodes, 'eclass1', eclass1, 'eclass2', eclass2, ...
    'observed', onodes);

prior = [sm.segmentationClusters.priorProbability];
prior = prior./sum(prior); %make sure prior probabilities sum to 1
bnet.CPD{1} = tabular_CPD(bnet, 1, 'CPT', prior);
mixmat = zeros([Q M]);
mu = zeros(O, Q, M);
sigma = zeros(O, O, Q, M);
[cm, ccv, cmr] = meanAndCovOfKnownPoints(sm.segmentationClusters, M);

for j = 1:length(sm.segmentationClusters)
    mixmat(j,:) = cmr{j};
    mu(:,j,:) = cm{j}';
    sigma(:,:,j,:) = ccv{j};
end

bnet.CPD{2} = tabular_CPD(bnet, 2, 'CPT', mixmat);
bnet.CPD{3} = gaussian_CPD(bnet, 3, 'mean', mu, 'cov', sigma);

transmat = zeros(Q,Q);

for i = 1:Q
    transmat(i,:) = prior .* sm.allowedTransitions(i,:) ./ (sum(prior.*sm.allowedTransitions(i,:)));
end
  
bnet.CPD{4} = tabular_CPD(bnet,4,'CPT', transmat, 'adjustable',0); %sm.allowedTransitions);
