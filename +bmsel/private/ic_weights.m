function w = ic_weights(ic)
%IC_WEIGHTS  Akaike (or Schwarz) weights from information-criterion values.
%
%   w = ic_weights(ic)     ic is a [1 x K] vector of summed AIC or BIC values.
%
%       Delta_k = IC_k - min_j IC_j
%       w_k     = exp(-Delta_k / 2) / sum_j exp(-Delta_j / 2)
%
%   w_k can be read as the probability that model k is the best of the set
%   (Burnham & Anderson, 2002). The weights sum to 1.

    d = ic - min(ic);
    w = exp(-0.5 * d);
    w = w / sum(w);
end
