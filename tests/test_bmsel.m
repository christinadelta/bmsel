function test_bmsel()
%TEST_BMSEL  Run all checks for the bmsel package.
%
%   Run from the repository root:   >> cd tests; test_bmsel
%   or from anywhere, if the repository root is on the path.
%
%   Tests 1-5 need only base MATLAB. Test 6 compares bmsel.bms with SPM's
%   spm_BMS and runs only if SPM is on the path.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));            % repository root, holds +bmsel

fprintf('\nbmsel tests\n-----------\n');

%% 1. Fixed regression values (computed independently in Python/SciPy)
% The deterministic parts of bms (alpha, exp_r, F1, F0, BOR) must match
% these reference numbers to high precision.
lme = [-100 -104 -110
        -95  -93  -99
       -120 -121 -118
        -80  -85  -90
       -101 -100 -105
        -99 -104 -103];
r = bmsel.bms(lme, 'Nsamp', 1e4, 'Seed', 1);
check(max(abs(r.alpha - [5.0306114726 2.3392789613 1.6301095660])) < 1e-6, ...
      'alpha matches reference');
check(max(abs(r.exp_r - [0.5589568303 0.2599198846 0.1811232851])) < 1e-8, ...
      'exp_r matches reference');
check(abs(r.F1 - (-596.2885758658)) < 1e-6, 'F1 matches reference');
check(abs(r.F0 - (-595.9248435304)) < 1e-6, 'F0 matches reference');
check(abs(r.bor - 0.5899436260) < 1e-8, 'BOR matches reference');

%% 2. Basic properties
check(abs(sum(r.exp_r) - 1) < 1e-12, 'exp_r sums to 1');
check(abs(sum(r.xp) - 1) < 1e-12,    'xp sums to 1');
check(abs(sum(r.pxp) - 1) < 1e-12,   'pxp sums to 1');
check(all(abs(sum(r.g, 2) - 1) < 1e-12), 'posterior assignments sum to 1 per subject');
check(r.bor >= 0 && r.bor <= 1, 'BOR is a probability');

%% 3. Identical evidence for all models -> nothing to choose
lme0 = repmat(-100, 30, 3);
r0 = bmsel.bms(lme0, 'Nsamp', 1e5, 'Seed', 1);
check(max(abs(r0.exp_r - 1/3)) < 1e-10, 'equal evidence gives equal frequencies');
check(r0.bor > 0.5, 'equal evidence gives high BOR');
check(max(abs(r0.pxp - 1/3)) < 0.02, 'equal evidence gives pxp near 1/K');

%% 4. Overwhelming evidence for model 1 -> model 1 wins outright
lme1 = repmat([0 -50 -50], 40, 1);
r1 = bmsel.bms(lme1, 'Nsamp', 1e5, 'Seed', 1);
check(r1.exp_r(1) > 0.9, 'strong evidence gives high frequency for the winner');
check(r1.xp(1) > 0.999, 'strong evidence gives xp near 1');
check(r1.bor < 1e-6, 'strong evidence gives BOR near 0');
check(r1.pxp(1) > 0.999, 'strong evidence gives pxp near 1');

%% 5. compare(): shapes, consistency and reproducibility
rng(7);
Nsub = 25; K = 3;
NLL = 100 + 10 * randn(Nsub, K);
k   = [2 3 4];
n   = randi([150 250], Nsub, 1);
c = bmsel.compare(NLL, k, n, 'ModelNames', {'A', 'B', 'C'}, ...
                  'Nsamp', 1e5, 'Seed', 3, 'Verbose', false);
check(isequal(size(c.AIC), [Nsub K]) && isequal(size(c.BIC), [Nsub K]), 'AIC/BIC have the right size');
check(max(max(abs(c.AIC - (2*NLL + 2*repmat(k, Nsub, 1))))) < 1e-10, 'AIC formula');
check(max(max(abs(c.BIC - (2*NLL + repmat(k, Nsub, 1) .* log(n))))) < 1e-10, 'BIC formula');
check(abs(sum(c.wAIC) - 1) < 1e-12 && abs(sum(c.wBIC) - 1) < 1e-12, 'IC weights sum to 1');
check(sum(c.cntAIC) == Nsub && sum(c.cntBIC) == Nsub, 'winner counts sum to Nsub');
check(isequal(c.lme, -0.5 * c.BIC), 'default evidence is -BIC/2');
check(height(c.tables.group) == K && height(c.tables.subject) == Nsub, 'tables have the right size');
c2 = bmsel.compare(NLL, k, n, 'ModelNames', {'A', 'B', 'C'}, ...
                   'Nsamp', 1e5, 'Seed', 3, 'Verbose', false);
check(isequal(c.xp, c2.xp), 'same seed gives the same xp');
ca = bmsel.compare(NLL, k, n, 'Evidence', 'aic', 'Nsamp', 1e4, 'Verbose', false);
check(isequal(ca.lme, -0.5 * ca.AIC), 'Evidence=aic uses -AIC/2');
cu = bmsel.compare(NLL, k, n, 'Evidence', -NLL, 'Nsamp', 1e4, 'Verbose', false);
check(isequal(cu.lme, -NLL), 'user evidence matrix is used as given');
check(isscalar(bmsel.compare(NLL, k, 200, 'Nsamp', 1e4, 'Verbose', false).bor), 'scalar n is accepted');

% a Dirichlet with a small alpha (< 1) exercises the a < 1 branch of the
% gamma sampler
xs = bmsel.bms([0 -1; -1 0; 0 0], 'Alpha0', [0.5 0.5], 'Nsamp', 1e4, 'Seed', 1);
check(abs(sum(xs.xp) - 1) < 1e-12, 'alpha < 1 sampling works');

%% 6. Against SPM (only if SPM is installed)
if exist('spm_BMS', 'file') == 2
    rng(11);
    L = -100 + 5 * randn(40, 4);
    L(:, 2) = L(:, 2) + 3;                       % make model 2 somewhat better
    [a_spm, er_spm, xp_spm, pxp_spm, bor_spm] = spm_BMS(L, 1e6, 0, 0, 1);
    rb = bmsel.bms(L, 'Nsamp', 1e6, 'Seed', 1);
    check(max(abs(rb.alpha - a_spm)) < 1e-3,  'alpha matches spm_BMS');
    check(max(abs(rb.exp_r - er_spm)) < 1e-4, 'exp_r matches spm_BMS');
    check(abs(rb.bor - bor_spm) < 1e-3,       'BOR matches spm_BMS');
    check(max(abs(rb.xp - xp_spm)) < 0.01,    'xp matches spm_BMS (Monte-Carlo tolerance)');
    check(max(abs(rb.pxp - pxp_spm)) < 0.01,  'pxp matches spm_BMS (Monte-Carlo tolerance)');
else
    fprintf('  (skipped) SPM not on the path: spm_BMS comparison not run\n');
end

fprintf('\nAll tests passed.\n\n');
end

function check(cond, msg)
    if cond
        fprintf('  ok    %s\n', msg);
    else
        error('bmsel:test:failed', 'FAILED: %s', msg);
    end
end
