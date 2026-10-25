function g = gamrnd_mt(a, n, stream)
%GAMRND_MT  Gamma(shape a, scale 1) random numbers without the Statistics Toolbox.
%
%   g = gamrnd_mt(a, n, stream)   returns an [n x 1] vector.
%
%   Uses the rejection method of Marsaglia & Tsang (2000). For shape a >= 1:
%
%       d = a - 1/3,   c = 1 / sqrt(9 d)
%       draw x ~ Normal(0,1),  v = (1 + c x)^3
%       draw u ~ Uniform(0,1)
%       accept  d * v  if  v > 0  and  ln u < x^2/2 + d - d v + d ln v
%
%   For shape a < 1 we use the boost
%
%       Gamma(a) = Gamma(a + 1) * U^(1/a),   U ~ Uniform(0,1)
%
%   Reference: Marsaglia, G., Tsang, W.W. (2000). A simple method for
%   generating gamma variables. ACM Transactions on Mathematical Software,
%   26(3), 363-372.

    if a < 1
        g = gamrnd_mt(a + 1, n, stream) .* (rand(stream, n, 1) .^ (1 / a));
        return;
    end
    d = a - 1/3;
    c = 1 / sqrt(9 * d);
    g = zeros(n, 1);
    need = true(n, 1);
    while any(need)
        idx = find(need);
        m   = numel(idx);
        x   = randn(stream, m, 1);
        v   = (1 + c * x) .^ 3;
        u   = rand(stream, m, 1);
        acc = (v > 0) & (log(u) < 0.5 * x.^2 + d - d * v + d * log(max(v, realmin)));
        g(idx(acc))    = d * v(acc);
        need(idx(acc)) = false;
    end
end
