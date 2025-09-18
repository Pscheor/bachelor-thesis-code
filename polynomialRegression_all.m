% === Initialization ===
VPs = unique(trial.VPCode);
numBins = 6;
latVar = 'Pupil_base_avg';

% Store individual regression parameters
intercepts_lin = [];
betas_lin = [];

intercepts_quad = [];
betas_quad1 = [];
betas_quad2 = [];

% Store all (x, y) points
x_all = [];
y_all = [];

for i = 1:numel(VPs)
    vpData = trial(strcmp(trial.VPCode, VPs(i)), :);
    pupil_data = vpData.(latVar);
    codes = vpData.SDTCode;

    if numel(pupil_data) < numBins || all(isnan(pupil_data))
        continue
    end

    % Binning
    edges = [-inf, quantile(pupil_data, (1:numBins-1)/numBins), inf];
    binIdx = discretize(pupil_data, edges);

    dprimes = NaN(1, numBins);
    pupilMeans = NaN(1, numBins);

    for b = 1:numBins
        idx = binIdx == b;
        if any(idx)
            pupilMeans(b) = mean(pupil_data(idx), 'omitnan');
            [dprimes(b), ~] = computeSDT(codes(idx));
        end
    end

    if all(~isnan(pupilMeans)) && all(~isnan(dprimes))
        % Center pupil values per subject
        pupilMeans_demeaned = pupilMeans - mean(pupilMeans, 'omitnan');

        % Linear regression
        mdl_lin = fitlm(pupilMeans_demeaned, dprimes, 'linear');
        intercepts_lin(end+1) = mdl_lin.Coefficients.Estimate(1);
        betas_lin(end+1) = mdl_lin.Coefficients.Estimate(2);

        % Quadratic regression
        mdl_quad = fitlm(pupilMeans_demeaned, dprimes, 'quadratic');
        intercepts_quad(end+1) = mdl_quad.Coefficients.Estimate(1);
        betas_quad1(end+1) = mdl_quad.Coefficients.Estimate(2);
        betas_quad2(end+1) = mdl_quad.Coefficients.Estimate(3);

        % Store individual data points
        x_all = [x_all; pupilMeans_demeaned(:)];
        y_all = [y_all; dprimes(:)];
    end
end

% Remove NaNs
valid = ~isnan(x_all) & ~isnan(y_all);
x = x_all(valid);
y = y_all(valid);

% === Average parameters across participants ===
m_b0_lin  = mean(intercepts_lin, 'omitnan');
m_b1_lin  = mean(betas_lin, 'omitnan');

m_b0_quad = mean(intercepts_quad, 'omitnan');
m_b1_quad = mean(betas_quad1, 'omitnan');
m_b2_quad = mean(betas_quad2, 'omitnan');

% === Compute regression lines from averaged parameters ===
xplot = linspace(min(x), max(x), 100);
y_lin = m_b0_lin + m_b1_lin * xplot;
y_quad = m_b0_quad + m_b1_quad * xplot + m_b2_quad * xplot.^2;

% === Model evaluation (based on real data) ===
n = length(y);
SStot = sum((y - mean(y)).^2);

% Linear
y_pred_lin = m_b0_lin + m_b1_lin * x;
SSres_lin = sum((y - y_pred_lin).^2);
R2_lin = 1 - SSres_lin / SStot;
k_lin = 2;
AIC_lin = n * log(SSres_lin / n) + 2 * k_lin;
BIC_lin = n * log(SSres_lin / n) + k_lin * log(n);

% Quadratic
y_pred_quad = m_b0_quad + m_b1_quad * x + m_b2_quad * x.^2;
SSres_quad = sum((y - y_pred_quad).^2);
R2_quad = 1 - SSres_quad / SStot;
k_quad = 3;
AIC_quad = n * log(SSres_quad / n) + 2 * k_quad;
BIC_quad = n * log(SSres_quad / n) + k_quad * log(n);

% === Plot ===
figure; hold on;

% Original data points
plot(x, y, 'k.', 'MarkerSize', 6);

% Regression curves (dashed)
plot(xplot, y_lin, 'b--', 'LineWidth', 2);
plot(xplot, y_quad, 'r--', 'LineWidth', 2);

xlabel('Centered pupil mean');
ylabel('Performance d′');
title('Regression using averaged parameters (across participants)');

% Legend with one-line stats per model
legend({...
    'Data points', ...
    sprintf('Linear (--): R² = %.3f, AIC = %.1f, BIC = %.1f', R2_lin, AIC_lin, BIC_lin), ...
    sprintf('Quadratic (--): R² = %.3f, AIC = %.1f, BIC = %.1f', R2_quad, AIC_quad, BIC_quad) ...
    }, ...
    'Location', 'best', ...
    'FontSize', 9, ...
    'Interpreter', 'none');

grid on;
