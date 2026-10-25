function xp = dirichlet_exceedance(alpha, Nsamp, stream)
%DIRICHLET_EXCEEDANCE  Exceedance probabilities of a Dirichlet distribution.
%
%   xp = dirichlet_exceedance(alpha, Nsamp, stream)
%
%   xp(k) = P( r_k > r_j for all j ~= k ),   r ~ Dirichlet(alpha)
%
%   Computed by Monte-Carlo sampling. A Dirichlet(alpha) sample is obtained
%   from independent Gamma draws:
%
%       z_k ~ Gamma(shape = alpha_k, scale = 1),   r_k = z_k / sum_j z_j
%
%   Because the normalisation is the same for every k, the model with the
%   largest z_k is also the one with the largest r_k, so we do not need to
%   divide. xp(k) is then the fraction of samples in which k has the maximum.
%
%   Samples are drawn in blocks so that memory stays small for large Nsamp.
%
%   stream : a RandStream (base MATLAB) for reproducibility.

    K      = numel(alpha);
    block  = 2e5;
    counts = zeros(1, K);
    done   = 0;
    while done < Nsamp
        n = min(block, Nsamp - done);
        z = zeros(n, K);
        for k = 1:K
            z(:, k) = gamrnd_mt(alpha(k), n, stream);
        end
        [~, w] = max(z, [], 2);
        counts = counts + accumarray(w, 1, [K 1])';
        done = done + n;
    end
    xp = counts / Nsamp;
end
