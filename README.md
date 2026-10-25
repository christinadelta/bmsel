# bmsel

**Bayesian model selection for group studies, in one call. Base MATLAB only.**

`bmsel` takes the fit results you already have (one negative log-likelihood per
subject and model, the number of free parameters, and the number of trials) and
returns everything you normally report when you compare computational models
across a group of participants:

- AIC and BIC per subject and summed over subjects
- Akaike and Schwarz weights
- how many subjects each model fits best
- random-effects Bayesian model selection (Stephan et al., 2009): expected
  model frequencies, exceedance probabilities, **protected** exceedance
  probabilities and the Bayesian omnibus risk (Rigoux et al., 2014)
- two result tables (group and per subject), a figure, and CSV/MAT output if
  you want it

It needs no toolboxes: no SPM, no VBA, no Statistics Toolbox. The random-effects
part is checked against SPM's `spm_BMS` (see [Tests](#tests)).

## Why another tool?

The same statistics are available in [SPM](https://www.fil.ion.ucl.ac.uk/spm/)
(`spm_BMS`), the [VBA toolbox](https://mbb-team.github.io/VBA-toolbox/)
(`VBA_groupBMC`) and [cbm](https://github.com/payampiray/cbm). Those are large
packages, or tied to their own fitting procedures. `bmsel` is small, readable,
dependency-free and does the whole fixed-effects and random-effects comparison in
one function. Every calculation is documented with its maths, both in the code
comments and in the [Methods](#methods) section below, so you can see exactly
what is being computed.

## Installation

Download or clone the repository and add its root folder (the one that contains
`+bmsel`) to your MATLAB path:

```matlab
addpath('path/to/bmsel');
```

Requires MATLAB R2016b or later.

## Quick start

```matlab
% NLL : [subjects x models] negative log-likelihood at the best-fitting parameters
% k   : [1 x models] number of free parameters
% n   : trials per subject (scalar, or one value per subject)

res = bmsel.compare(NLL, [2 3 3], n, 'ModelNames', {'M1', 'M2', 'M3'});
bmsel.plot(res);
```

This prints two tables and returns a struct `res` with every quantity (see
below). To also save CSV tables, a figure and a `.mat` file:

```matlab
res = bmsel.compare(NLL, k, n, 'ModelNames', names, 'OutDir', 'results', 'Seed', 1);
```

If you already have log model evidence from somewhere else (a Laplace
approximation, a variational free energy, a cross-validated log-likelihood),
skip the information criteria and run the random-effects part directly:

```matlab
rfx = bmsel.bms(lme);        % lme is [subjects x models]
```

Run `examples/example_simulated.m` for a worked example. Note that it fits
nothing: it makes up an NLL matrix for three pretend models, so you can see the
input format, the output, and that the built-in model frequencies are recovered.

## Functions

| Function | What it does |
|---|---|
| `bmsel.compare(NLL, k, n, ...)` | Full comparison: AIC, BIC, weights, winner counts, random-effects BMS, tables, optional files |
| `bmsel.bms(lme, ...)` | Random-effects BMS only, from a log-evidence matrix |
| `bmsel.plot(res, ...)` | Bar plots of expected frequencies and protected exceedance probabilities |

Type `help bmsel.compare` (and so on) in MATLAB for all options.

### Options of `bmsel.compare`

| Option | Default | Meaning |
|---|---|---|
| `'ModelNames'` | `{'M1','M2',...}` | names used in tables and plots |
| `'Evidence'` | `'bic'` | log evidence for the random-effects part: `'bic'`, `'aic'`, or your own `[subjects x models]` matrix |
| `'Nsamp'` | `1e6` | Dirichlet samples for the exceedance probabilities |
| `'Seed'` | `[]` | random seed, for reproducible exceedance probabilities |
| `'OutDir'` | `''` | if set, writes `<prefix>_group_summary.csv`, `<prefix>_per_subject.csv`, `<prefix>_results.mat` and `<prefix>_bms.png` there |
| `'Prefix'` | `'bmsel'` | file-name prefix for the saved files |
| `'Verbose'` | `true` | print the tables |

### Output of `bmsel.compare`

| Field | Size | Meaning |
|---|---|---|
| `AIC`, `BIC` | subjects × models | information criteria per subject |
| `sumAIC`, `sumBIC`, `dAIC`, `dBIC` | 1 × models | summed criteria, and difference from the best model |
| `wAIC`, `wBIC` | 1 × models | Akaike and Schwarz weights |
| `cntAIC`, `cntBIC` | 1 × models | number of subjects each model fits best |
| `lme` | subjects × models | the log evidence passed to the random-effects analysis |
| `alpha` | 1 × models | posterior Dirichlet counts |
| `exp_r` | 1 × models | expected model frequencies |
| `xp` | 1 × models | exceedance probabilities |
| `pxp` | 1 × models | protected exceedance probabilities (**the number to report**) |
| `bor` | scalar | Bayesian omnibus risk |
| `g` | subjects × models | posterior probability that each subject used each model |
| `tables.group`, `tables.subject` | tables | everything above, ready for `writetable` |

### How to read the results

- **Fixed effects** (summed AIC/BIC, weights) assume that everyone in the group
  uses the same model. They are easy to interpret but are driven by a few
  subjects with large differences.
- **Random effects** (`exp_r`, `xp`, `pxp`, `bor`) allow different people to use
  different models and ask how frequent each model is in the population. This is
  the recommended analysis for group studies.
- `exp_r` is the estimated proportion of people best described by each model.
- `xp` is the probability that a model is more frequent than all the others.
  It goes to 1 for the leading model as the sample grows, even if the models
  are nearly equally frequent, so it can overstate the evidence.
- `bor` is the probability that all models are equally frequent (the "null").
  If it is high (above about 0.25) the data do not really separate the models.
- `pxp` is `xp` corrected for that risk, and is the quantity to report.

## Methods

Notation: $`i = 1, \dots, N`$ indexes subjects and $`k = 1, \dots, K`$ indexes
models. $`y_i`$ is subject $`i`$'s data, $`\mathrm{NLL}_{ik}`$ the negative
log-likelihood of $`y_i`$ under model $`k`$ at the best-fitting parameters,
$`p_k`$ the number of free parameters of model $`k`$, and $`n_i`$ the number of data
points (trials) of subject $`i`$. All logarithms are natural.

### 1. Information criteria

$$
\mathrm{AIC}_{ik} = 2\,\mathrm{NLL}_{ik} + 2\,p_k
\qquad
\mathrm{BIC}_{ik} = 2\,\mathrm{NLL}_{ik} + p_k \ln n_i
$$

(Akaike, 1974; Schwarz, 1978). Lower is better. Both are on the deviance
scale, so a difference of about 10 between two models is usually taken as
strong evidence. BIC penalises parameters more than AIC once $`n_i > e^2 \approx 7`$.

### 2. Fixed-effects summaries

The group criterion is the sum over subjects, $`\mathrm{IC}_k = \sum_i \mathrm{IC}_{ik}`$,
and $`\Delta_k = \mathrm{IC}_k - \min_j \mathrm{IC}_j`$. The Akaike (or Schwarz)
weight of model $`k`$ is

$$
w_k = \frac{\exp(-\Delta_k / 2)}{\sum_j \exp(-\Delta_j / 2)}
$$

(Burnham & Anderson, 2002). $`w_k`$ is the probability that model $`k`$ is the
best of the set, under the assumption that all subjects use the same model.
`bmsel` also counts, for each model, the number of subjects for which it has the
lowest AIC and the lowest BIC.

### 3. Log model evidence

The random-effects analysis needs the log marginal likelihood (log model
evidence) $`\ln p(y_i \mid k)`$ of each subject under each model. The default
uses the BIC approximation (Kass & Raftery, 1995):

$$
\ln p(y_i \mid k) \approx -\frac{1}{2}\mathrm{BIC}_{ik}
= -\mathrm{NLL}_{ik} - \frac{p_k}{2} \ln n_i
$$

With `'Evidence', 'aic'` it uses $`-\frac{1}{2}\,\mathrm{AIC}_{ik}`$ instead, and
with a numeric `'Evidence'` matrix you can supply any other estimate (Laplace
approximation, variational free energy, cross-validated log-likelihood). Call
the resulting matrix $`L_{ik} = \ln p(y_i \mid k)`$.

### 4. Random-effects model: the generative story

Following Stephan et al. (2009), each subject uses one model, $`m_i \in \{1,\dots,K\}`$,
drawn from a population with unknown model frequencies
$`r = (r_1, \dots, r_K)`$, $`\sum_k r_k = 1`$:

$$
r \sim \mathrm{Dirichlet}(\alpha_0), \qquad
m_i \mid r \sim \mathrm{Categorical}(r), \qquad
p(y_i \mid m_i = k) = \exp(L_{ik})
$$

The prior is flat, $`\alpha_0 = (1, \dots, 1)`$. We want the posterior over $`r`$.

### 5. Variational inference

The posterior is approximated by $`q(r)\,\prod_i q(m_i)`$ with
$`q(r) = \mathrm{Dirichlet}(\alpha)`$ and $`q(m_i = k) = g_{ik}`$. The optimal
factors are found by iterating until $`\alpha`$ stops changing:

$$
\mathbb{E}_q[\ln r_k] = \psi(\alpha_k) - \psi\left(\sum_j \alpha_j\right)
$$

$$
g_{ik} = \frac{\exp(L_{ik} + \mathbb{E}_q[\ln r_k])}
              {\sum_j \exp(L_{ij} + \mathbb{E}_q[\ln r_j])}
$$

$$
\alpha_k = \alpha_{0k} + \sum_{i=1}^{N} g_{ik}
$$

where $`\psi`$ is the digamma function. $`g_{ik}`$ is the posterior probability
that subject $`i`$'s data were generated by model $`k`$ (returned as `g`), and
$`\alpha_k`$ counts the "effective" number of subjects using model $`k`$ plus
the prior count.

### 6. Expected frequencies and exceedance probabilities

The expected model frequency is the Dirichlet mean,

$$
\mathbb{E}[r_k \mid y] = \frac{\alpha_k}{\sum_j \alpha_j}
$$

(returned as `exp_r`). The exceedance probability of model $`k`$ is the
posterior probability that it is more frequent than every other model,

$$
\varphi_k = P(r_k > r_j \ \text{for all } j \neq k \mid y)
$$

(returned as `xp`). There is no closed form for $`K > 2`$, so it is estimated by
Monte-Carlo sampling: draw $`r \sim \mathrm{Dirichlet}(\alpha)`$ many times
(default $`10^6`$) and count how often each $`r_k`$ is the largest. A Dirichlet
sample is built from independent Gamma draws,
$`z_k \sim \mathrm{Gamma}(\alpha_k, 1)`$, $`r_k = z_k / \sum_j z_j`$; the Gamma
samples use the method of Marsaglia & Tsang (2000), so no Statistics Toolbox
is needed.

### 7. Bayesian omnibus risk and protected exceedance probabilities

Exceedance probabilities can be misleading: with enough subjects the leading
model's $`\varphi_k`$ approaches 1 even when all models are almost equally
frequent. Rigoux et al. (2014) therefore compare the model above ($`H_1`$:
frequencies differ) against the null $`H_0`$: $`r_k = 1/K`$ for all $`k`$.

Under $`H_0`$ each subject's model is a fair draw, so the log evidence of
$`H_0`$ is

$$
F_0 = \sum_{i=1}^{N} \ln \left( \frac{1}{K} \sum_{k=1}^{K} \exp(L_{ik}) \right)
$$

The log evidence of $`H_1`$ is approximated by the variational free energy at
the solution of step 5 (the same quantity as in SPM's `spm_BMS_F`):

$$
F_1 = A + B - C
$$

where

$$
A = \sum_{i=1}^{N}\sum_{k=1}^{K} g_{ik}\,(L_{ik} + \mathbb{E}_q[\ln r_k]) - \sum_{i=1}^{N}\sum_{k=1}^{K} g_{ik} \ln g_{ik}
$$

$$
B = \ln\Gamma\left(\sum_k \alpha_{0k}\right) - \sum_k \ln\Gamma(\alpha_{0k}) + \sum_k (\alpha_{0k} - 1)\,\mathbb{E}_q[\ln r_k]
$$

$$
C = \ln\Gamma\left(\sum_k \alpha_k\right) - \sum_k \ln\Gamma(\alpha_k) + \sum_k (\alpha_k - 1)\,\mathbb{E}_q[\ln r_k]
$$

$`A`$ is the expected log joint of data and model labels plus the entropy of
$`q(m)`$; $`B`$ is $`\mathbb{E}_q[\ln p(r)]`$ and $`C`$ is $`\mathbb{E}_q[\ln q(r)]`$.

With equal prior probability on the two hypotheses, the **Bayesian omnibus
risk** is the posterior probability of the null,

$$
\mathrm{BOR} = p(H_0 \mid y) = \frac{1}{1 + \exp(F_1 - F_0)}
$$

and the **protected exceedance probability** averages $`\varphi_k`$ over the
two hypotheses (under $`H_0`$ every model is "best" with probability $`1/K`$):

$$
\tilde{\varphi}_k = (1 - \mathrm{BOR})\,\varphi_k + \frac{\mathrm{BOR}}{K}
$$

(returned as `pxp`). This is the quantity to report.

## Tests

```matlab
cd tests
test_bmsel
```

The tests check the deterministic quantities ($`\alpha`$, `exp_r`, $`F_0`$,
$`F_1`$, BOR) against reference values computed independently in Python/SciPy,
check the limiting cases (equal evidence, overwhelming evidence), and check the
formulas, sizes and reproducibility of `bmsel.compare`. If SPM is on the path,
`bmsel.bms` is also compared with `spm_BMS` on the same log-evidence matrix.

## Citing

If you use `bmsel`, please cite the repository (see `CITATION.cff`; a DOI is
provided through Zenodo) and the methods papers:

- Stephan, K. E., Penny, W. D., Daunizeau, J., Moran, R. J., & Friston, K. J.
  (2009). Bayesian model selection for group studies. *NeuroImage, 46*(4),
  1004–1017.
- Rigoux, L., Stephan, K. E., Friston, K. J., & Daunizeau, J. (2014). Bayesian
  model selection for group studies — revisited. *NeuroImage, 84*, 971–985.

## References

- Akaike, H. (1974). A new look at the statistical model identification.
  *IEEE Transactions on Automatic Control, 19*(6), 716–723.
- Burnham, K. P., & Anderson, D. R. (2002). *Model selection and multimodel
  inference* (2nd ed.). Springer.
- Kass, R. E., & Raftery, A. E. (1995). Bayes factors. *Journal of the
  American Statistical Association, 90*(430), 773–795.
- Marsaglia, G., & Tsang, W. W. (2000). A simple method for generating gamma
  variables. *ACM Transactions on Mathematical Software, 26*(3), 363–372.
- Rigoux, L., Stephan, K. E., Friston, K. J., & Daunizeau, J. (2014). Bayesian
  model selection for group studies — revisited. *NeuroImage, 84*, 971–985.
- Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of
  Statistics, 6*(2), 461–464.
- Stephan, K. E., Penny, W. D., Daunizeau, J., Moran, R. J., & Friston, K. J.
  (2009). Bayesian model selection for group studies. *NeuroImage, 46*(4),
  1004–1017.

## Licence

MIT. See `LICENSE`.
