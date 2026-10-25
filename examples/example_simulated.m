%% bmsel example: synthetic fit results for three models
%
% IMPORTANT: nothing is fitted in this example. There are no real models and
% no task data. The script makes up an NLL matrix directly (one negative
% log-likelihood per subject and model), so that you can see:
%   - the input format that bmsel.compare expects (NLL, k, n)
%   - what the output looks like
%   - that the estimated model frequencies recover the ones we built in
%
% In a real analysis, NLL comes from fitting each model to each subject's
% data (for example with fminsearch, fmincon or bayesopt), k is the number
% of free parameters of each model, and n is the number of trials per
% subject.
%
% Run from the repository root, or add the repository root to the path.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));          % repository root, holds +bmsel

rng(2026);

%% Make up fit results
% 60 subjects and three models, M1, M2, M3. We decide that in the population
% 60% of people are best described by M1, 30% by M2 and 10% by M3. For each
% subject, the "true" model gets the lowest NLL; the other two are worse by a
% random amount of 3 to 9 log-likelihood units.
Nsub   = 60;
names  = {'M1', 'M2', 'M3'};
k      = [2 3 3];                       % pretend numbers of free parameters
n      = randi([180 240], Nsub, 1);     % pretend trials per subject
r_true = [0.6 0.3 0.1];                 % true model frequencies

true_model = sum(rand(Nsub, 1) > cumsum(r_true), 2) + 1;   % one label per subject
NLL = zeros(Nsub, 3);
for i = 1:Nsub
    base = 0.55 * n(i) + 8 * randn;                % NLL of the true model
    penalty = 3 + 6 * rand(1, 3);                  % how much worse the others are
    penalty(true_model(i)) = 0;
    NLL(i, :) = base + penalty;
end

%% Compare the models
res = bmsel.compare(NLL, k, n, 'ModelNames', names, 'Seed', 1);

%% Plot
bmsel.plot(res, 'Title', 'synthetic example');

%% Look at the tables
disp(res.tables.group);
disp(head(res.tables.subject));

%% Did we recover the frequencies we built in?
fprintf('True frequencies:      %s\n', mat2str(r_true, 2));
fprintf('Estimated frequencies: %s\n', mat2str(res.exp_r, 2));

%% Only the random-effects part, from your own log evidence
% If you already have log model evidence (for example from a Laplace
% approximation or a cross-validated log-likelihood), skip the information
% criteria and call bms directly:
lme = -0.5 * res.BIC;                   % here we just reuse -BIC/2
rfx = bmsel.bms(lme, 'Seed', 1);
fprintf('pxp from bms alone:    %s\n', mat2str(rfx.pxp, 3));
