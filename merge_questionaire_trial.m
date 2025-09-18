%% === Daten laden
load('C:\Users\larsd\Desktop\Bachelorarbeit\PupilPersonality\behavior_trial.mat');
load('C:\Users\larsd\Desktop\Bachelorarbeit\PupilPersonality\questionaireData_all.mat');
questionaireData_all = summaryTable;

%% === VP-Codes vereinheitlichen
behavior_trial.VPCode = upper(string(behavior_trial.VPCode));
questionaireData_all.VP_Code = upper(string(questionaireData_all.VP_Code));

%% === Nur gewünschte Fragebogen-Spalten auswählen
varsToKeep = { ...
    'VP_Code', ...
    'Agree_Mean', 'Consc_Mean', 'Extra_Mean', ...
    'Neuro_Mean', 'Open_Mean', ...
    'OCI_Sum', 'Autism_Sum', ...
    'Age', 'Education', 'Gender', 'Hearing' ...
};

questionaireData_reduced = questionaireData_all(:, varsToKeep);

%% === Mergen über VPCode
% Rename VP_Code -> VPCode zur Übereinstimmung
questionaireData_reduced.Properties.VariableNames{1} = 'VPCode';

% Join durchführen
behavior_trial = outerjoin(behavior_trial, questionaireData_reduced, ...
    'Keys', 'VPCode', ...
    'MergeKeys', true, ...
    'Type', 'left');

%% === Speichern
trial_BQ = behavior_trial;
save('C:\Users\larsd\Desktop\Final\trial_BQ.mat', 'trial_BQ');
fprintf('trial_BQ gespeichert.\n');

% === Fehlende Fragebogendaten prüfen
% Eine Spalte aus dem Fragebogenbereich nehmen, z. B. "Agree_Mean"
missingQuestionnaireIdx = isnan(behavior_trial.Agree_Mean);

% Fehlende VPCode(s) extrahieren
missingVPs = unique(behavior_trial.VPCode(missingQuestionnaireIdx));

% Ausgabe (nur wenn welche fehlen)
if ~isempty(missingVPs)
    fprintf('\n⚠️  Keine Fragebogendaten gefunden für folgende VPCode(s):\n');
    disp(missingVPs);
else
    fprintf('\n✅ Alle VPCodes haben passende Fragebogendaten.\n');
end
