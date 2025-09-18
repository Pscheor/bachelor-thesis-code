%% === Laden der Daten
load('C:\Users\larsd\Desktop\Final\trial_pupil.mat');   % enthält trial_pupil
load('C:\Users\larsd\Desktop\Final\trial_BQ.mat');      % enthält trial_BQ

%% === Sicherstellen, dass Trialnmb numerisch ist
trial_pupil.Trialnmb = double(trial_pupil.Trialnmb);
trial_BQ.Trialnmb = double(trial_BQ.TrialNumber);
%% === Schritt 1: Block berechnen in trial_pupil
blockSize = 50;
trial_pupil.BlockNumber = ceil(trial_pupil.Trialnmb / blockSize);

%% === Schritt 2: Entferne erste 2 Trials je Block
toRemove = false(height(trial_pupil), 1);
vpList = unique(trial_pupil.VPCode);

for i = 1:numel(vpList)
    vp = vpList{i};
    for cond = ["liberal", "conservative"]
        % Alle vorhandenen Blöcke
        idx_vp_cond = strcmp(trial_pupil.VPCode, vp) & strcmpi(trial_pupil.Cond, cond);
        blocks = unique(trial_pupil.BlockNumber(idx_vp_cond));

        for b = blocks'
            idx_block = idx_vp_cond & trial_pupil.BlockNumber == b;

            trialNums = trial_pupil.Trialnmb(idx_block);
            [~, sortIdx] = sort(trialNums);
            idx_block_rows = find(idx_block);

            % max 2 Trials entfernen
            toRemove(idx_block_rows(sortIdx(1:min(2,end)))) = true;
        end
    end
end

% Bereinigte Tabelle
trial_pupil_clean = trial_pupil(~toRemove, :);

%% === Schritt 3: Neue TrialNumber pro Block generieren
trial_pupil_clean.TrialNumber = NaN(height(trial_pupil_clean), 1);

for i = 1:numel(vpList)
    vp = vpList{i};
    for cond = ["liberal", "conservative"]
        idx_vp_cond = strcmp(trial_pupil_clean.VPCode, vp) & strcmpi(trial_pupil_clean.Cond, cond);
        blocks = unique(trial_pupil_clean.BlockNumber(idx_vp_cond));

        for b = blocks'
            idx_block = idx_vp_cond & trial_pupil_clean.BlockNumber == b;
            idx_block_rows = find(idx_block);
            trial_pupil_clean.TrialNumber(idx_block_rows) = (1:length(idx_block_rows))';
        end
    end
end

%% === Schritt 4: Condition-Kleinschreibung angleichen
trial_pupil_clean.Condition = regexprep(lower(trial_pupil_clean.Cond), '^(.)', '${upper($1)}');

% Du kannst auch die alte Spalte überschreiben:
trial_pupil_clean.Cond = [];

%% === Schritt 5: Join-Vorbereitung
% Spaltennamen angleichen
trial_pupil_clean.Properties.VariableNames{'Condition'} = 'Condition';
trial_pupil_clean.Properties.VariableNames{'TrialNumber'} = 'TrialNumber';

% Sicherstellen, dass Join-Schlüssel vorhanden sind
requiredKeys = {'VPCode', 'Condition', 'TrialNumber'};
assert(all(ismember(requiredKeys, trial_pupil_clean.Properties.VariableNames)), 'Fehlende Join-Spalten in trial_pupil_clean');
assert(all(ismember(requiredKeys, trial_BQ.Properties.VariableNames)), 'Fehlende Join-Spalten in trial_BQ');

%% === Join durchführen
% Vor dem Merge: sicherstellen, dass trial_pupil_clean keine Duplikate hat
[~, uniqueIdx] = unique(trial_pupil_clean(:, {'VPCode', 'Condition', 'TrialNumber'}), 'rows');
trial_pupil_clean = trial_pupil_clean(uniqueIdx, :);

% Merge mit eindeutiger Zuordnung
trial = innerjoin(trial_BQ, trial_pupil_clean, ...
    'Keys', {'VPCode', 'Condition', 'TrialNumber'});


%% === Prüfung: Fehlende Zuordnungen
missingPupil = isnan(trial.Pupil_base_left);
n_missing = sum(missingPupil);

if n_missing > 0
    fprintf('\n⚠️  %d Zeilen konnten keine Pupillendaten zugeordnet werden.\n', n_missing);
    disp(unique(trial.VPCode(missingPupil)));
else
    fprintf('\n✅ Alle Trials wurden korrekt gemerged.\n');
end

%%sortieren
% Definiere Reihenfolge für 'Condition': Liberal vor Conservative
trial.Condition = categorical(trial.Condition, {'Liberal', 'Conservative'}, 'Ordinal', true);

trial.BlockNumber = trial.BlockNumber_trial_BQ;
    trial.BlockNumber_trial_BQ = [];
    trial.BlockNumber_trial_pupil_clean = [];
    trial.Trialnmb_trial_BQ = [];
    trial.Trialnmb_trial_pupil_clean = [];
% Neue Spaltenreihenfolge:
% - VPCode, Condition, TrialNumber, BlockNumber, dann der Rest
varsFixed = {'VPCode', 'Condition', 'BlockNumber', 'TrialNumber'};

% Alle restlichen Spalten, ohne die ersten vier
varsRest = setdiff(trial.Properties.VariableNames, varsFixed, 'stable');

% Endgültige Spaltenreihenfolge
trial = trial(:, [varsFixed, varsRest]);


    % Sortiere nach VP, Condition, Block, TrialNumber
trial = sortrows(trial, {'VPCode', 'Condition', 'BlockNumber', 'TrialNumber'});

%% === Speichern
save('C:\Users\larsd\Desktop\Final\trial.mat', 'trial');
fprintf('\n📄 Merged-Tabelle gespeichert als: trial.mat\n');
