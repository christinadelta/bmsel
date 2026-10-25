function h = plot(res, varargin)
%PLOT  Bar plots of the random-effects model selection results.
%
%   h = bmsel.plot(res)
%   h = bmsel.plot(res, 'Name', Value, ...)
%
%   res is the output of bmsel.compare or bmsel.bms. Two panels are drawn:
%   the expected model frequencies (exp_r) and the protected exceedance
%   probabilities (pxp), with the BOR in the title.
%
%   OPTIONS
%     'ModelNames'  labels for the x axis (default: res.models if present)
%     'Title'       text added before the panel titles (default '')
%     'Visible'     'on' (default) or 'off'
%
%   Returns the figure handle.

p = inputParser;
p.FunctionName = 'bmsel.plot';
addRequired(p, 'res', @isstruct);
addParameter(p, 'ModelNames', {}, @(x) iscellstr(x) || isstring(x));
addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Visible', 'on', @(x) any(strcmpi(x, {'on', 'off'})));
parse(p, res, varargin{:});
opt = p.Results;

K = numel(res.exp_r);
if ~isempty(opt.ModelNames)
    names = cellstr(opt.ModelNames);
elseif isfield(res, 'models')
    names = res.models;
else
    names = arrayfun(@(i) sprintf('M%d', i), 1:K, 'UniformOutput', false);
end
ttl = char(opt.Title);
if ~isempty(ttl), ttl = [ttl ': ']; end

h = figure('Visible', opt.Visible, 'Color', 'w', 'Position', [100 100 760 320]);

subplot(1, 2, 1);
bar(res.exp_r, 'FaceColor', [0.30 0.45 0.70]);
ylim([0 1]); xlim([0.5 K + 0.5]);
set(gca, 'XTick', 1:K, 'XTickLabel', names);
ylabel('expected frequency  E[r_k]');
title([ttl 'model frequencies']);
box off;

subplot(1, 2, 2);
bar(res.pxp, 'FaceColor', [0.70 0.35 0.30]);
ylim([0 1]); xlim([0.5 K + 0.5]);
set(gca, 'XTick', 1:K, 'XTickLabel', names);
ylabel('protected exceedance probability');
title(sprintf('%sPXP  (BOR = %.3f)', ttl, res.bor));
box off;
end
