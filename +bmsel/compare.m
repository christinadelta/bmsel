function res = compare(NLL, k, n, varargin)
%COMPARE  Group-level model comparison from per-subject fit results.
%
%   res = bmsel.compare(NLL, k, n)
%   res = bmsel.compare(NLL, k, n, 'Name', Value, ...)
%
%   One call gives you everything you normally report when comparing
%   computational models across subjects:
%     - AIC and BIC per subject, and summed over subjects
%     - Akaike and Schwarz weights
%     - how many subjects each model wins
%     - random-effects Bayesian model selection: expected model frequencies,
%       exceedance and protected exceedance probabilities, and the BOR
%   Only base MATLAB is needed.
%
%   INPUT
%     NLL   [Nsub x K] negative log-likelihood of each subject's data under
%           each model, at the best-fitting parameters (natural log).
%     k     [1 x K] number of free parameters of each model.
%     n     number of data points (trials, choices) per subject: either a
%           scalar, if it is the same for everyone, or an [Nsub x 1] vector.
%
%   OPTIONS (name-value pairs)
%     'ModelNames'  cell array of K names (default {'M1','M2',...})
%     'Evidence'    what to use as log model evidence for the random-effects
%                   analysis: 'bic' (default), 'aic', or an [Nsub x K]
%                   matrix of your own log evidence (Laplace, variational,
%                   cross-validated, ...).
%     'Nsamp'       Dirichlet samples for exceedance probabilities (1e6)
%     'Seed'        random seed for reproducible xp (default [])
%     'OutDir'      folder to write results to. If given, two CSV tables, a
%                   figure (PNG) and a .mat file are saved there (default '').
%     'Prefix'      file-name prefix for the saved files (default 'bmsel')
%     'Verbose'     print the result tables to the Command Window (true)
%
%   OUTPUT  res, a struct with fields
%     models, k, n, Nsub, K          the inputs
%     NLL, AIC, BIC                  [Nsub x K]
%     sumNLL, sumAIC, sumBIC         [1 x K] summed over subjects
%     dAIC, dBIC                     [1 x K] difference from the best model
%     wAIC, wBIC                     [1 x K] Akaike and Schwarz weights
%     bestAIC, bestBIC               [Nsub x 1] index of each subject's best model
%     cntAIC, cntBIC                 [1 x K] number of subjects won
%     evidence, lme                  what was used for the random-effects part
%     alpha, exp_r, xp, pxp, bor, g  random-effects results (see bmsel.bms)
%     tables.group, tables.subject   the two result tables
%
%   EXAMPLE
%     NLL = [120 118 125; 95 99 94; 110 109 115];   % 3 subjects, 3 models
%     res = bmsel.compare(NLL, [2 3 3], 200, 'ModelNames', {'RW','RW+','HGF'});
%     bmsel.plot(res);
%
%   MATHS (see README.md for the full write-up)
%     AIC_ik = 2 NLL_ik + 2 k_k
%     BIC_ik = 2 NLL_ik + k_k ln n_i
%     weights:   w_k = exp(-Delta_k/2) / sum_j exp(-Delta_j/2)
%     evidence:  ln p(y_i | k) ~ -BIC_ik / 2   (or -AIC_ik / 2)
%     then bmsel.bms on the evidence matrix.
%
%   See also bmsel.bms, bmsel.plot.

% ---------------------------------------------------------------- inputs --
p = inputParser;
p.FunctionName = 'bmsel.compare';
addRequired(p, 'NLL', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
addRequired(p, 'k',   @(x) isnumeric(x) && isvector(x) && all(x >= 0));
addRequired(p, 'n',   @(x) isnumeric(x) && all(x(:) > 0));
addParameter(p, 'ModelNames', {},    @(x) iscellstr(x) || isstring(x));
addParameter(p, 'Evidence',   'bic', @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'Nsamp',      1e6,   @(x) isscalar(x) && x >= 1);
addParameter(p, 'Seed',       [],    @(x) isempty(x) || (isscalar(x) && x >= 0));
addParameter(p, 'OutDir',     '',    @(x) ischar(x) || isstring(x));
addParameter(p, 'Prefix',     'bmsel', @(x) ischar(x) || isstring(x));
addParameter(p, 'Verbose',    true,  @(x) islogical(x) || isnumeric(x));
parse(p, NLL, k, n, varargin{:});
opt = p.Results;

NLL = double(NLL);
[Nsub, K] = size(NLL);
k = double(k(:))';
if numel(k) ~= K
    error('bmsel:compare:kSize', 'k must have one entry per model (column of NLL).');
end
if isscalar(n)
    n = repmat(double(n), Nsub, 1);
else
    n = double(n(:));
    if numel(n) ~= Nsub
        error('bmsel:compare:nSize', 'n must be a scalar or have one entry per subject (row of NLL).');
    end
end
if any(~isfinite(NLL(:)))
    error('bmsel:compare:nonFinite', ...
        'NLL contains NaN or Inf. Remove those subjects first (all models must have a fit for every subject).');
end
if K < 2
    error('bmsel:compare:oneModel', 'Need at least two models to compare.');
end
if isempty(opt.ModelNames)
    names = arrayfun(@(i) sprintf('M%d', i), 1:K, 'UniformOutput', false);
else
    names = cellstr(opt.ModelNames); names = names(:)';
    if numel(names) ~= K
        error('bmsel:compare:namesSize', 'ModelNames must have one name per model.');
    end
end
outdir = char(opt.OutDir);
prefix = char(opt.Prefix);

% ---------------------------------------------------------- AIC and BIC --
% For subject i and model k, with NLL the negative log-likelihood at the
% best-fitting parameters, k_k free parameters and n_i data points:
%
%   AIC_ik = 2 NLL_ik + 2 k_k                 (Akaike, 1974)
%   BIC_ik = 2 NLL_ik + k_k ln(n_i)           (Schwarz, 1978)
%
% Lower is better. BIC punishes extra parameters more than AIC once
% n_i > e^2 (about 7 trials), so BIC tends to prefer simpler models.
Kmat = repmat(k, Nsub, 1);
AIC  = 2 * NLL + 2 * Kmat;
BIC  = 2 * NLL + Kmat .* log(n);          % log(n) is Nsub x 1, expands over models

% -------------------------------------------- fixed-effects group summary --
% Summing over subjects treats everyone as using the same model (a
% "fixed-effects" view). Differences in summed IC are on the deviance
% scale: a difference of 10 or more is usually taken as strong.
sumNLL = sum(NLL, 1);
sumAIC = sum(AIC, 1);
sumBIC = sum(BIC, 1);
dAIC   = sumAIC - min(sumAIC);
dBIC   = sumBIC - min(sumBIC);

% Akaike / Schwarz weights (Burnham & Anderson, 2002):
%   w_k = exp(-Delta_k / 2) / sum_j exp(-Delta_j / 2),  Delta_k = IC_k - min IC
% w_k is the (fixed-effects) probability that model k is the best of the set.
wAIC = ic_weights(sumAIC);
wBIC = ic_weights(sumBIC);

% --------------------------------------------------- per-subject winners --
% A simple count of how many subjects each model fits best. It ignores how
% big the differences are, so it is a descriptive number only.
[~, bestAIC] = min(AIC, [], 2);
[~, bestBIC] = min(BIC, [], 2);
cntAIC = accumarray(bestAIC, 1, [K 1])';
cntBIC = accumarray(bestBIC, 1, [K 1])';

% --------------------------------------- random-effects model selection --
% The random-effects analysis needs the log model evidence ln p(y_i | k).
% With the BIC approximation (Schwarz, 1978; Kass & Raftery, 1995):
%
%   ln p(y_i | k) ~ -BIC_ik / 2  =  -NLL_ik - (k_k / 2) ln n_i
%
% The AIC version is -AIC_ik / 2. Both are rough. If you have a better
% estimate (Laplace approximation, variational free energy, cross-validated
% log-likelihood) pass it with 'Evidence'.
if isnumeric(opt.Evidence)
    lme = double(opt.Evidence);
    if ~isequal(size(lme), [Nsub K])
        error('bmsel:compare:evidenceSize', 'Evidence matrix must be the same size as NLL.');
    end
    evidence = 'user';
else
    switch lower(char(opt.Evidence))
        case 'bic', lme = -0.5 * BIC; evidence = 'bic';
        case 'aic', lme = -0.5 * AIC; evidence = 'aic';
        otherwise
            error('bmsel:compare:evidence', 'Evidence must be ''bic'', ''aic'' or a numeric matrix.');
    end
end
rfx = bmsel.bms(lme, 'Nsamp', opt.Nsamp, 'Seed', opt.Seed);

% ---------------------------------------------------------------- tables --
grpT = table((1:K)', names(:), k(:), sumNLL(:), sumAIC(:), sumBIC(:), ...
             dAIC(:), dBIC(:), wAIC(:), wBIC(:), cntAIC(:), cntBIC(:), ...
             rfx.alpha(:), rfx.exp_r(:), rfx.xp(:), rfx.pxp(:), repmat(rfx.bor, K, 1), ...
    'VariableNames', {'model_idx', 'model', 'k', 'sumNLL', 'sumAIC', 'sumBIC', ...
             'dAIC', 'dBIC', 'wAIC', 'wBIC', 'n_bestAIC', 'n_bestBIC', ...
             'alpha', 'exp_r', 'xp', 'pxp', 'BOR'});

vn   = matlab.lang.makeValidName(names);     % safe column names (no spaces etc.)
subT = table((1:Nsub)', n, 'VariableNames', {'subject', 'n'});
for m = 1:K, subT.(['NLL_' vn{m}]) = NLL(:, m); end
for m = 1:K, subT.(['AIC_' vn{m}]) = AIC(:, m); end
for m = 1:K, subT.(['BIC_' vn{m}]) = BIC(:, m); end
for m = 1:K, subT.(['p_' vn{m}])   = rfx.g(:, m); end   % posterior model probability per subject
subT.best_AIC = names(bestAIC)';
subT.best_BIC = names(bestBIC)';

% ---------------------------------------------------------------- output --
res = struct();
res.models   = names;
res.k        = k;
res.n        = n;
res.Nsub     = Nsub;
res.K        = K;
res.NLL      = NLL;
res.AIC      = AIC;
res.BIC      = BIC;
res.sumNLL   = sumNLL;
res.sumAIC   = sumAIC;
res.sumBIC   = sumBIC;
res.dAIC     = dAIC;
res.dBIC     = dBIC;
res.wAIC     = wAIC;
res.wBIC     = wBIC;
res.bestAIC  = bestAIC;
res.bestBIC  = bestBIC;
res.cntAIC   = cntAIC;
res.cntBIC   = cntBIC;
res.evidence = evidence;
res.lme      = lme;
res.alpha    = rfx.alpha;
res.exp_r    = rfx.exp_r;
res.xp       = rfx.xp;
res.pxp      = rfx.pxp;
res.bor      = rfx.bor;
res.g        = rfx.g;
res.F1       = rfx.F1;
res.F0       = rfx.F0;
res.tables   = struct('group', grpT, 'subject', subT);

% ----------------------------------------------------------------- print --
if opt.Verbose
    fprintf('\nbmsel: %d subjects, %d models, evidence = %s\n', Nsub, K, evidence);
    fprintf('\n  Fixed effects (lower IC = better):\n');
    fprintf('  %-12s %3s %10s %10s %10s %8s %8s %7s %7s\n', ...
        'Model', 'k', 'sumNLL', 'sumAIC', 'sumBIC', 'dAIC', 'dBIC', 'wAIC', 'wBIC');
    for m = 1:K
        fprintf('  %-12s %3d %10.1f %10.1f %10.1f %8.1f %8.1f %7.3f %7.3f\n', ...
            names{m}, k(m), sumNLL(m), sumAIC(m), sumBIC(m), dAIC(m), dBIC(m), wAIC(m), wBIC(m));
    end
    fprintf('\n  Random effects (BOR = %.3f):\n', rfx.bor);
    fprintf('  %-12s %9s %9s %8s %7s %7s %7s\n', ...
        'Model', '#bestAIC', '#bestBIC', 'alpha', 'exp_r', 'xp', 'pxp');
    for m = 1:K
        fprintf('  %-12s %9d %9d %8.2f %7.3f %7.3f %7.3f\n', ...
            names{m}, cntAIC(m), cntBIC(m), rfx.alpha(m), rfx.exp_r(m), rfx.xp(m), rfx.pxp(m));
    end
    fprintf('\n');
end

% ------------------------------------------------------------------ save --
if ~isempty(outdir)
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    writetable(grpT, fullfile(outdir, [prefix '_group_summary.csv']));
    writetable(subT, fullfile(outdir, [prefix '_per_subject.csv']));
    save(fullfile(outdir, [prefix '_results.mat']), 'res');
    try
        h = bmsel.plot(res, 'Visible', 'off');
        saveas(h, fullfile(outdir, [prefix '_bms.png']));
        close(h);
    catch ME
        warning('bmsel:compare:figure', 'Could not save the figure: %s', ME.message);
    end
    if opt.Verbose
        fprintf('  Results written to %s\n\n', outdir);
    end
end
end
