% %% Vorbereitung
% clear;
% dataPath = 'C:\Uni\BA\Daten\preproc';
% load(fullfile(dataPath, 'trial_cleaned.mat'));  % lädt 'trial'
% trial.VPCode = string(trial.VPCode);

% Nutzung der vorhandenen Spalte
latVar = 'Pupil_base_avg';  % Pupillenmaß
numBins = 6;

% Gruppen, Bedingungen
conds = unique(trial.Condition);
groups = ["introvert", "extravert"];

% Farben
colors = struct( ...
    'introvert', [0.2 0.4 1], ...
    'extravert', [0.9 0.2 0.2], ...
    'Liberal', [0.2 0.6 1], ...
    'Conservative', [0.9 0.2 0.2], ...
    'introvert_Liberal', [0.1 0.4 0.9], ...
    'introvert_Conservative', [0.1 0.2 0.6], ...
    'extravert_Liberal', [1 0.4 0.4], ...
    'extravert_Conservative', [0.7 0.1 0.1]);

%% 1. Gesamtauswertung über alle VPs
VPs = unique(trial.VPCode);
dMat_all = [];
xMat_all = [];

for i = 1:numel(VPs)
    vp = VPs(i);
    vpData = trial(trial.VPCode == vp, :);
    pupil_data = vpData.(latVar);
    codes = vpData.SDTCode;

    edges = [-inf, quantile(pupil_data, (1:numBins-1)/numBins), inf];
    binIdx = discretize(pupil_data, edges);

    dprimes = NaN(1, numBins);
    pupilMeans = NaN(1, numBins);

    for b = 1:numBins
        thisBin = binIdx == b;
        if any(thisBin)
            pupilMeans(b) = mean(pupil_data(thisBin), 'omitnan');
            [dprimes(b), ~] = computeSDT(codes(thisBin));
        end
    end

    dMat_all(end+1, :) = dprimes;
    xMat_all(end+1, :) = pupilMeans;
end

% Plot Gesamt
figure('Name', 'Gesamt d′ vs Pupillengröße', 'Units', 'normalized', 'Position', [0.3 0.3 0.45 0.5]);
meanX = mean(xMat_all, 1, 'omitnan');
meanY = mean(dMat_all, 1, 'omitnan');
SEM = std(dMat_all - mean(dMat_all, 2, 'omitnan'), 0, 1, 'omitnan') / sqrt(size(dMat_all,1));

errorbar(meanX, meanY, SEM, '-o', ...
    'LineWidth', 2, 'Color', [0.2 0.2 0.2], 'MarkerFaceColor', [0.5 0.5 0.5]);
xlabel('Pupil size (mean per bin)'); ylabel('Performance d′');
title('Performance d′ vs Pupilsize – across all');
grid on; ylim padded;

%% 2. Extraversion Split (ExtraversionGroup2)
extraData = []; extraX = [];
introData = []; introX = [];

validGroups = ["introvert", "extravert"];
trial = trial(ismember(trial.ExtraversionGroup2, validGroups), :);
VPs = unique(trial.VPCode);

for i = 1:numel(VPs)
    vp = VPs(i);
    data = trial(trial.VPCode == vp, :);
    pupil_data = data.(latVar);
    codes = data.SDTCode;
    group = unique(data.ExtraversionGroup2);
    if numel(group) ~= 1, continue; end

    edges = [-inf, quantile(pupil_data, (1:numBins-1)/numBins), inf];
    binIdx = discretize(pupil_data, edges);

    dprimes = NaN(1, numBins);
    pupilMeans = NaN(1, numBins);

    for b = 1:numBins
        thisBin = binIdx == b;
        if any(thisBin)
            pupilMeans(b) = mean(pupil_data(thisBin), 'omitnan');
            [dprimes(b), ~] = computeSDT(codes(thisBin));
        end
    end

    if group == "extravert"
        extraData(end+1, :) = dprimes;
        extraX(end+1, :) = pupilMeans;
    elseif group == "introvert"
        introData(end+1, :) = dprimes;
        introX(end+1, :) = pupilMeans;
    end
end

% Plot: Extraversion-Split
figure('Name', 'Performance d′ vs Pupil size - Extraversion split', 'Units', 'normalized', 'Position', [0.3 0.3 0.45 0.5]); hold on;
meanX = mean(extraX, 1, 'omitnan');
meanY = mean(extraData, 1, 'omitnan');
semY  = std(extraData - mean(extraData,2,'omitnan'),0,1,'omitnan') / sqrt(size(extraData,1));
errorbar(meanX, meanY, semY, '-o', 'Color', colors.extravert, 'MarkerFaceColor', colors.extravert*0.8, 'LineWidth',2, 'DisplayName','Extravert');

meanX = mean(introX, 1, 'omitnan');
meanY = mean(introData, 1, 'omitnan');
semY  = std(introData - mean(introData,2,'omitnan'),0,1,'omitnan') / sqrt(size(introData,1));
errorbar(meanX, meanY, semY, '-o', 'Color', colors.introvert, 'MarkerFaceColor', colors.introvert*0.8, 'LineWidth',2, 'DisplayName','Introvert');

xlabel('Pupil size (mean per bin)'); ylabel('Performance d′');
title('Performance d′ vs Pupil size – Extraversion split');
legend('Location','best'); grid on; ylim padded;

%% 3. Condition Split
allResults = struct();
for c = 1:numel(conds)
    cond = conds(c);
    data = trial(trial.Condition == cond, :);
    VPs = unique(data.VPCode);
    dMat = []; xMat = [];

    for i = 1:numel(VPs)
        vp = VPs(i);
        vpData = data(data.VPCode == vp, :);
        pupil_data = vpData.(latVar);
        codes = vpData.SDTCode;

        edges = [-inf, quantile(pupil_data, (1:numBins-1)/numBins), inf];
        binIdx = discretize(pupil_data, edges);

        dprimes = NaN(1, numBins);
        pupilMeans = NaN(1, numBins);

        for b = 1:numBins
            thisBin = binIdx == b;
            if any(thisBin)
                pupilMeans(b) = mean(pupil_data(thisBin), 'omitnan');
                [dprimes(b), ~] = computeSDT(codes(thisBin));
            end
        end

        dMat(end+1,:) = dprimes;
        xMat(end+1,:) = pupilMeans;
    end

    allResults.(char(cond)).dMat = dMat;
    allResults.(char(cond)).xMat = xMat;
end

% Plot
figure('Name','d′ vs Pupillengröße – nach Condition','Units','normalized','Position',[0.3 0.3 0.45 0.5]); hold on;
for c = 1:numel(conds)
    cond = conds(c);
    dMat = allResults.(char(cond)).dMat;
    xMat = allResults.(char(cond)).xMat;

    meanX = mean(xMat, 1, 'omitnan');
    meanY = mean(dMat, 1, 'omitnan');
    semY = std(dMat - mean(dMat,2,'omitnan'), 0, 1, 'omitnan') / sqrt(size(dMat,1));

    errorbar(meanX, meanY, semY, '-o', 'Color', colors.(char(cond)), ...
        'MarkerFaceColor', colors.(char(cond))*0.8, 'LineWidth',2, 'DisplayName', char(cond));
end
xlabel('Pupil size (mean per bin)'); ylabel('Performance d′');
title('Performance d′ vs Pupil size – Condition spllit');
legend('Location','best'); grid on; ylim padded;

%% 4. Extraversion × Condition Split
% Struktur für Ergebnisse
allRes = struct();
for g = 1:numel(groups)
    group = groups(g);
    gData = trial(trial.ExtraversionGroup2 == group, :);

    for c = 1:numel(conds)
        cond = conds(c);
        data = gData(gData.Condition == cond, :);
        VPs = unique(data.VPCode);

        dMat = []; xMat = [];

        for i = 1:numel(VPs)
            vp = VPs(i);
            vpData = data(data.VPCode == vp, :);
            pupil_data = vpData.(latVar);
            codes = vpData.SDTCode;

            edges = [-inf, quantile(pupil_data, (1:numBins-1)/numBins), inf];
            binIdx = discretize(pupil_data, edges);

            dprimes = NaN(1, numBins);
            pupilMeans = NaN(1, numBins);

            for b = 1:numBins
                thisBin = binIdx == b;
                if any(thisBin)
                    pupilMeans(b) = mean(pupil_data(thisBin), 'omitnan');
                    [dprimes(b), ~] = computeSDT(codes(thisBin));
                end
            end

            dMat(end+1,:) = dprimes;
            xMat(end+1,:) = pupilMeans;
        end

        allRes.(group).(char(cond)).dMat = dMat;
        allRes.(group).(char(cond)).xMat = xMat;
    end
end

% Plot
figure('Name','d′ vs Pupillengröße – Extraversion × Condition', ...
       'Units','normalized','Position',[0.3 0.3 0.55 0.6]); hold on;

for g = 1:numel(groups)
    group = groups(g);
    for c = 1:numel(conds)
        cond = conds(c);
        dMat = allRes.(group).(char(cond)).dMat;
        xMat = allRes.(group).(char(cond)).xMat;

        meanX = mean(xMat, 1, 'omitnan');
        meanY = mean(dMat, 1, 'omitnan');
        semY = std(dMat - mean(dMat,2,'omitnan'), 0, 1, 'omitnan') / sqrt(size(dMat,1));

        color = colors.([char(group) '_' char(cond)]);
        errorbar(meanX, meanY, semY, '-o', 'Color', color, ...
            'MarkerFaceColor', color*0.8, 'LineWidth',2, ...
            'DisplayName', sprintf('%s – %s', group, cond));
    end
end
xlabel('Pupil size (mean per bin)'); ylabel('Performance d′');
title('Performance d′ vs Pupil size – Extraversion & Condition');
legend('Location','best'); grid on; ylim padded;
