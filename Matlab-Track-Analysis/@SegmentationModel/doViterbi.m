function guessedStates = doViterbi (sm, track)
%function guessedStates = doViterbi (sm, track)

[data, known_points] = sm.createHMMData(track);
B = mixgauss_prob(data{1}, sm.hmm_mean, sm.hmm_cov, sm.hmm_mix);

inds = find(known_points{1} > 0);
B(:,inds) = eps;
B(sub2ind(size(B), known_points{1}(inds), inds)) = 1;

%where data(:,:,ex) is OxT where O is the size of the observation vector. Finally, use

[guessedStates] = viterbi_path(sm.hmm_prior, sm.hmm_transmat, B);
