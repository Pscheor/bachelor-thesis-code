% Ensure VPCode is string
trial.VPCode = string(trial.VPCode);

% Compute mean extraversion score per participant
vp_stats = varfun(@mean, trial, 'InputVariables', 'Extra_Mean', ...
    'GroupingVariables', {'VPCode', 'ExtraversionGroup2'});

% Extract groups
groups = unique(vp_stats.ExtraversionGroup2);

% Print descriptive statistics
fprintf('\n=== Descriptive Statistics for ExtraversionGroup2 ===\n');
for i = 1:numel(groups)
    group = groups(i);
    data = vp_stats.mean_Extra_Mean(vp_stats.ExtraversionGroup2 == group);

    n = numel(data);
    m = mean(data, 'omitnan');
    s = std(data, 'omitnan');
    med = median(data, 'omitnan');
    minVal = min(data, [], 'omitnan');
    maxVal = max(data, [], 'omitnan');

    fprintf('%s:\n', group);
    fprintf('  n      = %d\n', n);
    fprintf('  Mean   = %.2f\n', m);
    fprintf('  SD     = %.2f\n', s);
    fprintf('  Median = %.2f\n', med);
    fprintf('  Min    = %.2f\n', minVal);
    fprintf('  Max    = %.2f\n\n', maxVal);
end

% Liste aller VPs
VPs = unique(trial.VPCode);

% Initialisierung
c_values = NaN(numel(VPs),1);

for i = 1:numel(VPs)
    vp = VPs(i);
    vp_data = trial(trial.VPCode == vp, :);
    sdt_codes = vp_data.SDTCode;

    % Berechne d′ und c mit deiner computeSDT Funktion
    [~, c] = computeSDT(sdt_codes);
    c_values(i) = c;
end

% Deskriptive Statistik für c
n      = numel(c_values);
M      = mean(c_values, 'omitnan');
SD     = std(c_values, 'omitnan');
med    = median(c_values, 'omitnan');
minVal = min(c_values, [], 'omitnan');
maxVal = max(c_values, [], 'omitnan');

% Ausgabe
fprintf('\n=== Descriptive Statistics for Response Criterion c (per VP, all trials) ===\n');
fprintf('n      = %d\n', n);
fprintf('Mean   = %.3f\n', M);
fprintf('SD     = %.3f\n', SD);
fprintf('Median = %.3f\n', med);
fprintf('Min    = %.3f\n', minVal);
fprintf('Max    = %.3f\n', maxVal);

% Liste aller VPs
VPs = unique(trial.VPCode);

% Leere Liste für d′-Werte
dprimes = NaN(numel(VPs), 1);

% d′-Berechnung pro VP über alle Trials
for i = 1:numel(VPs)
    vp = VPs(i);
    sdt = trial.SDTCode(trial.VPCode == vp);

    [d, ~] = computeSDT(sdt);
    dprimes(i) = d;
end

% Deskriptivstatistik über alle VPs
n      = sum(~isnan(dprimes));
m      = mean(dprimes, 'omitnan');
s      = std(dprimes, 'omitnan');
med    = median(dprimes, 'omitnan');
minVal = min(dprimes, [], 'omitnan');
maxVal = max(dprimes, [], 'omitnan');

% Ausgabe
fprintf('\n=== Descriptive Statistics for d′ (per VP, all trials) ===\n');
fprintf('n      = %d\n', n);
fprintf('Mean   = %.3f\n', m);
fprintf('SD     = %.3f\n', s);
fprintf('Median = %.3f\n', med);
fprintf('Min    = %.3f\n', minVal);
fprintf('Max    = %.3f\n', maxVal);

% Mittelwert der Pupil_base_avg pro VP berechnen
vp_stats = varfun(@mean, trial, 'InputVariables', 'Pupil_base_avg', ...
                  'GroupingVariables', 'VPCode');

% Die Spalte mit den Mittelwerten extrahieren
data = vp_stats.mean_Pupil_base_avg;

% Deskriptivstatistiken berechnen
n      = numel(data);                      % Anzahl VPs
m      = mean(data, 'omitnan');
s      = std(data, 'omitnan');
med    = median(data, 'omitnan');
minVal = min(data, [], 'omitnan');
maxVal = max(data, [], 'omitnan');

% Ausgabe
fprintf('\n=== Descriptive Statistics for Pupil_base_avg (averaged per VP) ===\n');
fprintf('n      = %d\n', n);
fprintf('Mean   = %.3f\n', m);
fprintf('SD     = %.3f\n', s);
fprintf('Median = %.3f\n', med);
fprintf('Min    = %.3f\n', minVal);
fprintf('Max    = %.3f\n', maxVal);

% Stelle sicher, dass VPCode als string vorliegt
trial.VPCode = string(trial.VPCode);

% Mittelwert der ReactionTime pro VP berechnen
vp_stats = varfun(@mean, trial, 'InputVariables', 'ReactionTime', ...
                  'GroupingVariables', 'VPCode');

% Extrahiere die mittleren RTs pro VP
RT_values = vp_stats.mean_ReactionTime;

% Deskriptive Statistik
n      = numel(RT_values);
M      = mean(RT_values, 'omitnan');
SD     = std(RT_values, 'omitnan');
med    = median(RT_values, 'omitnan');
minVal = min(RT_values, [], 'omitnan');
maxVal = max(RT_values, [], 'omitnan');

% Ausgabe
fprintf('\n=== Descriptive Statistics for Reaction Time (per VP) ===\n');
fprintf('n      = %d\n', n);
fprintf('Mean   = %.3f ms\n', M);
fprintf('SD     = %.3f ms\n', SD);
fprintf('Median = %.3f ms\n', med);
fprintf('Min    = %.3f ms\n', minVal);
fprintf('Max    = %.3f ms\n', maxVal);
 % Stelle sicher, dass VPCode als string vorliegt
trial.VPCode = string(trial.VPCode);

% Alter pro VP extrahieren (ein Wert pro Person)
vp_stats = varfun(@mean, trial, 'InputVariables', 'Age', ...
                  'GroupingVariables', 'VPCode');

% Alterswerte
ages = vp_stats.mean_Age;

% Deskriptive Statistik
n      = numel(ages);
M      = mean(ages, 'omitnan');
SD     = std(ages, 'omitnan');
med    = median(ages, 'omitnan');
minVal = min(ages, [], 'omitnan');
maxVal = max(ages, [], 'omitnan');

% Ausgabe
fprintf('\n=== Descriptive Statistics for Age (per VP) ===\n');
fprintf('n      = %d\n', n);
fprintf('Mean   = %.2f years\n', M);
fprintf('SD     = %.2f years\n', SD);
fprintf('Median = %.2f years\n', med);
fprintf('Min    = %.2f years\n', minVal);
fprintf('Max    = %.2f years\n', maxVal);
