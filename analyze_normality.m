function analyze_normality(trial)
    % Analyze normality of Extra_Mean per VP from trial table
    % Histogram, Q-Q plot, Shapiro-Wilk test (if available).

    % === Check variable ===
    if ~ismember('Extra_Mean', trial.Properties.VariableNames)
        error('Variable "Extra_Mean" not found in trial table.');
    end
    if ~ismember('VPCode', trial.Properties.VariableNames)
        error('Variable "VPCode" not found in trial table.');
    end

    % === Aggregate: one value per VP (mean) ===
    [uniqueVPs, ~, idx] = unique(trial.VPCode);
    data = splitapply(@mean, trial.Extra_Mean, idx);  % mean per VP
    data = data(~isnan(data)); % remove NaN

    % === Figure setup ===
    figure('Name','Normality Analysis: Extra_Mean', ...
           'Units','normalized','Position',[0.1 0.1 0.8 0.6]);
    tl = tiledlayout(1,2, 'TileSpacing','compact','Padding','compact');
    title(tl, 'Normality Check for Extra\_Mean (per VP)');
% === Histogram (bars touching) ===
ax1 = nexttile(tl,1);
[counts, edges] = histcounts(data, 'BinMethod','sturges');
centers = edges(1:end-1) + diff(edges)/2;

bar(ax1, centers, counts, 1, 'FaceColor',[0.2 0.5 0.8], 'EdgeColor','k'); 
xlabel(ax1, 'Extraversion score');   % <-- customize
ylabel(ax1, 'Number of participants');    % <-- customize
title(ax1, 'Histogram');
grid(ax1,'on'); box(ax1,'on');

% Auto axis limits
xlim(ax1, [min(data) max(data)] + range(data)*[-0.05 0.05]);
ylim(ax1, [0, max(counts)+1]);

    % === Q-Q Plot ===
    ax2 = nexttile(tl,2);
    qqplot(ax2, data);
    xlabel(ax2, 'Standard Normal Quantiles'); % <-- customize
    ylabel(ax2, 'Quantiles of Sample');      % <-- customize
    title(ax2, 'Q-Q Plot');
    grid(ax2,'on'); box(ax2,'on');

    % === Normality Tests ===
    fprintf('\n--- Normality Tests for Extra_Mean (per VP) ---\n');

    % Shapiro-Wilk (preferred)
    if exist('swtest','file')
        [hS,pS] = swtest(data);
        fprintf('Shapiro-Wilk test: H=%d, p=%.4f\n', hS, pS);
    else
        warning('Shapiro-Wilk test (swtest.m) not found. Please install from File Exchange.');
    end

    % Optional backup tests
    [hL,pL] = lillietest(data);
    fprintf('Lilliefors test: H=%d, p=%.4f\n', hL, pL);

    [hJ,pJ] = jbtest(data);
    fprintf('Jarque-Bera test: H=%d, p=%.4f\n', hJ, pJ);
end
