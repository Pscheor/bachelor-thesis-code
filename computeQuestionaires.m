% Pfade setzen
cd("C:\Users\larsd\Desktop\Final");
folder = ("C:\Users\larsd\Desktop\Final");
filename = 'data_PupilPersonality_2025-07-04_09-49';
filepath = fullfile(folder, filename);

% Speicherort festlegen
folderPath = "C:\Users\larsd\Desktop\Final"; % Anpassen!
fileName = 'questionaireData_all.mat';
fullPath = fullfile(folderPath, fileName);

% CSV-Datei einlesen
data = readtable('data_PupilPersonality_2025-07-04_09-49.csv');

% VP-Code & Basisinfos extrahieren
VP_Code   = data.ID01_01;
Age       = data.SD01_01;
Gender    = data.SD02; % bereits als numerisch angenommen
Education = data.SD03;
Hearing   = data.SD04_01;

% Skalenfelder definieren
Neuro_fields = {'P001_01','P001_02','P001_03','P001_04','P001_05','P001_06'};
Extra_fields = {'P002_01','P002_02','P002_03','P002_04','P002_05','P002_06'};
Open_fields  = {'P003_01','P003_02','P003_03','P003_04','P003_05','P003_06'};
Agree_fields = {'P004_01','P004_02','P004_03','P004_04','P004_05','P004_06'};
Cons_fields  = {'P005_01','P005_02','P005_03','P005_04','P005_05','P005_06'};
% Kombinierte Liste aller NEO-Items
NEO_fields = [Neuro_fields, Extra_fields, Open_fields, Agree_fields, Cons_fields];

OCI_fields   = [strcat('O001_', string(10:18)), strcat('O002_', compose('%02d', 1:9))];
Autism_fields = strcat('AQ01_', compose('%02d', 1:10));

% GDMS: Block 1-5 jeweils 5 Items
GDMS_fields = [strcat('GD01_', compose('%02d', 1:5)), ...
               strcat('GD02_', compose('%02d', 1:5)), ...
               strcat('GD03_', compose('%02d', 1:5)), ...
               strcat('GD04_', compose('%02d', 1:5)), ...
               strcat('GD05_', compose('%02d', 1:5))];
% === GDMS Subskalen definieren ===
GDMS_Rational     = {'GD01_01', 'GD02_01', 'GD03_01', 'GD04_01', 'GD05_01'};
GDMS_Intuitive    = {'GD01_02', 'GD02_02', 'GD03_02', 'GD04_02', 'GD05_02'};
GDMS_Dependent    = {'GD01_03', 'GD02_03', 'GD03_03', 'GD04_03', 'GD05_03'};
GDMS_Avoiding     = {'GD01_04', 'GD02_04', 'GD03_04', 'GD04_04', 'GD05_04'};
GDMS_Spontaneous  = {'GD01_05', 'GD02_05', 'GD03_05', 'GD04_05', 'GD05_05'};

% === Subskalenmittelwerte berechnen ===
GDMS_Rational_Mean    = mean(data{:, GDMS_Rational},    2, 'omitnan');
GDMS_Intuitive_Mean   = mean(data{:, GDMS_Intuitive},   2, 'omitnan');
GDMS_Dependent_Mean   = mean(data{:, GDMS_Dependent},   2, 'omitnan');
GDMS_Avoiding_Mean    = mean(data{:, GDMS_Avoiding},    2, 'omitnan');
GDMS_Spontaneous_Mean = mean(data{:, GDMS_Spontaneous}, 2, 'omitnan');

%% Rekodierung negativer Items (Skala: 1–5)
% Big5
% Umkodierung zu offiziellen Skala 0-4
for i = 1:length(NEO_fields)
    var = NEO_fields{i};
    if ismember(var, data.Properties.VariableNames)
        if max(data.(var), [], 'omitnan') > 4  % Schutz gegen doppelte Umkodierung
            data.(var) = data.(var) - 1;
        end
    end
end
% === Liste der negativ gepolten Items (NEO-FFI-30) ===
negItems = {'P002_04', ...
            'P003_01', 'P003_03', 'P003_05', ...
            'P004_01', 'P004_02', 'P004_03', 'P004_04', 'P004_06', ...
            'P005_06'};

% === Sichere Rekodierung (einmalig, verhindert doppelte Umkehrung)
for i = 1:length(negItems)
    var = negItems{i};
    if ismember(var, data.Properties.VariableNames)
        maxVal = max(data.(var), [], 'omitnan');
        minVal = min(data.(var), [], 'omitnan');

        % Nur rekodieren, wenn Daten im erwarteten Bereich (noch nicht rekodiert)
        if maxVal <= 4 && minVal >= 0 && maxVal >= 2
            data.(var) = 4 - data.(var);
            fprintf('Item %s wurde rekodiert.\n', var);
        else
            fprintf('Item %s wurde NICHT rekodiert (vermutlich bereits umkodiert).\n', var);
        end
    else
        warning('NEO-Item %s nicht gefunden und konnte nicht geprüft werden.', var);
    end
end

%% AQ-short Rekodierung negativ formulierter Items (Likert 1–4)
% Von 1–4 auf 0–3 umkodieren (vor Rekodierung negativer Items!)
Autism_fields = strcat('AQ01_', compose('%02d', 1:10));

for i = 1:length(Autism_fields)
    var = Autism_fields{i};
    if ismember(var, data.Properties.VariableNames)
        if max(data.(var), [], 'omitnan') >= 4  % Schutz gegen doppelte Umkodierung
            data.(var) = data.(var) - 1;
        end
    end
end

% Skala: 1 = stimme voll zu → 4 = stimme überhaupt nicht zu
% Umkodierung: 3 - Antwort
aqItemsToReverse = {'AQ01_02', 'AQ01_03', 'AQ01_04', 'AQ01_05', 'AQ01_06', 'AQ01_09'};

for i = 1:length(aqItemsToReverse)
    var = aqItemsToReverse{i};
    if ismember(var, data.Properties.VariableNames)
        data.(var) = 3 - data.(var);
    else
        warning('Item %s nicht in der Tabelle gefunden.', var);
    end
end

% ===AQ-Kurz Summenwert berechnen (Likert)
Autism_fields = strcat('AQ01_', compose('%02d', 1:10));
Autism_Sum= sum(data{:, Autism_fields}, 2, 'omitnan');  % alternativ: mean(...) für Mittelwert

% === Klassifikation basierend auf Likert-Summenwert
% Schwelle: >= 30 = auffällig
AutismBinary = Autism_Sum >= 30;

% Textliche Klassifikation
AutismClass = repmat("unauffällig", height(data), 1);
AutismClass(AutismBinary) = "auffällig";

% === Histogramm der AQ-Kurz-Summenwerte ===
figure;
histogram(Autism_Sum, 'BinWidth', 1, 'FaceColor', [0.2 0.6 1.0], 'EdgeColor', 'black');
xlabel('AQ-Kurz Likert-Summenwert (10–40)');
ylabel('Anzahl der Probanden');
title('Verteilung des AQ-Kurz-Summenwerts');

% Klassifikationslinie (Schwelle = 30)
hold on;
xline(30, '--r', 'Cut-off: 30', 'LineWidth', 2, 'LabelOrientation', 'horizontal');
grid on;

% Plot speichern (Ordner: 'plots', falls nicht vorhanden → erstellen)
plotFolder = fullfile(folder, 'plots');
if ~exist(plotFolder, 'dir')
    mkdir(plotFolder);
end

plotFile = fullfile(plotFolder, 'AQ_Kurz_Histogramm.png');
exportgraphics(gcf, plotFile, 'Resolution', 300);
close(gcf);  % Fenster schließen (für automatisierte Skripte)

disp(['AQ-Histogramm gespeichert unter: ', plotFile]);

%% OCI-Items von Skala 1–5 → 0–4 umwandeln 
% OCI-Items laut SoSci-Codierung
OCI_fields = [strcat('O001_', string(10:18)), strcat('O002_', compose('%02d', 1:9))];

% Umkodierung durchführen (nur falls Werte noch bei 1–5 liegen)
for i = 1:length(OCI_fields)
    var = OCI_fields{i};
    if ismember(var, data.Properties.VariableNames)
        % Nur umkodieren, wenn Werte > 0 (sonst schon angepasst)
        if max(data.(var), [], 'omitnan') > 4
            data.(var) = data.(var) - 1;
        end
    else
        warning('OCI-Item %s nicht gefunden.', var);
    end
end
OCI_Sum = sum(data{:, OCI_fields}, 2, 'omitnan');  % 0–72
OCI_Binary = OCI_Sum >= 21;
OCI_Class = repmat("unauffällig", height(data), 1);
OCI_Class(OCI_Binary) = "auffällig";

% === Histogramm für OCI (Gesamtwert) ===
figure;
histogram(OCI_Sum, 'BinWidth', 2, 'FaceColor', [0.4 0.6 0.8], 'EdgeColor', 'black');
xlabel('OCI-Gesamtwert (0–72)');
ylabel('Anzahl der Probanden');
title('Verteilung des OCI-Gesamtwerts');
grid on;

% Schwelle markieren (Cut-off für klinisch auffällig: 21 Punkte)
hold on;
xline(21, '--r', 'Cut-off: 21', 'LineWidth', 2, ...
    'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal');

% === Speichern im "plots" Ordner ===
plotFolder = fullfile(folder, 'plots');
if ~exist(plotFolder, 'dir')
    mkdir(plotFolder);
end

plotFile = fullfile(plotFolder, 'OCI_Gesamtwert_Histogramm.png');
exportgraphics(gcf, plotFile, 'Resolution', 300);
close(gcf);  % Fenster schließen, falls Skript automatisch läuft

disp(['OCI-Histogramm gespeichert unter: ', plotFile]);

%% Skalenmittelwerte berechnen
Neuro_Mean = mean(data{:, Neuro_fields}, 2, 'omitnan');
Extra_Mean = mean(data{:, Extra_fields}, 2, 'omitnan');
Open_Mean  = mean(data{:, Open_fields},  2, 'omitnan');
Agree_Mean = mean(data{:, Agree_fields}, 2, 'omitnan');
Consc_Mean = mean(data{:, Cons_fields},  2, 'omitnan');
GDMS_Mean  = mean(data{:, GDMS_fields},  2, 'omitnan');

% === Ausschluss: Probanden mit negativen Werten bei Big Five
negMask = Neuro_Mean < 0 | Extra_Mean < 0 | Open_Mean < 0 | ...
          Agree_Mean < 0 | Consc_Mean < 0;

% Anzahl merken und ausschließen
n_excluded = sum(negMask);
disp([num2str(n_excluded), ' Probanden wurden ausgeschlossen wegen negativer Big Five Werte.']);

% Alle betroffenen Variablen bereinigen
VP_Code    = VP_Code(~negMask);
Age        = Age(~negMask);
Gender     = Gender(~negMask);
Education  = Education(~negMask);
Hearing    = Hearing(~negMask);

Neuro_Mean = Neuro_Mean(~negMask);
Extra_Mean = Extra_Mean(~negMask);
Open_Mean  = Open_Mean(~negMask);
Agree_Mean = Agree_Mean(~negMask);
Consc_Mean = Consc_Mean(~negMask);

OCI_Sum    = OCI_Sum(~negMask);
OCI_Binary = OCI_Binary(~negMask);
OCI_Class  = OCI_Class(~negMask);

Autism_Sum     = Autism_Sum(~negMask);
AutismBinary = AutismBinary(~negMask);
AutismClass  = AutismClass(~negMask);

GDMS_Mean            = GDMS_Mean(~negMask);
GDMS_Rational_Mean   = GDMS_Rational_Mean(~negMask);
GDMS_Intuitive_Mean  = GDMS_Intuitive_Mean(~negMask);
GDMS_Dependent_Mean  = GDMS_Dependent_Mean(~negMask);
GDMS_Avoiding_Mean   = GDMS_Avoiding_Mean(~negMask);
GDMS_Spontaneous_Mean = GDMS_Spontaneous_Mean(~negMask);

% Zusammengefasste Tabelle: In summaryTable integrieren
summaryTable = table(VP_Code, Age, Gender, Education, Hearing, ...
    Neuro_Mean, Extra_Mean, Open_Mean, Agree_Mean, Consc_Mean, ...
    OCI_Sum, OCI_Binary, OCI_Class, ...
    Autism_Sum, AutismBinary, AutismClass, ...
    GDMS_Mean, ...
    GDMS_Rational_Mean, ...
    GDMS_Intuitive_Mean, ...
    GDMS_Dependent_Mean, ...
    GDMS_Avoiding_Mean, ...
    GDMS_Spontaneous_Mean);

% --- Spaltennamen aus der Tabelle entnehmen
varNames = summaryTable.Properties.VariableNames;

% --- Mittelwerte & Standardabweichungen berechnen (nur Big Five Skalen)
big5Vars = {'Neuro_Mean', 'Extra_Mean', 'Open_Mean', 'Agree_Mean', 'Consc_Mean'};

% --- Vektor mit NaN-Werten für alle anderen Variablen
nCols = width(summaryTable);
meanRow = nan(1, nCols);
stdRow  = nan(1, nCols);

% --- Fülle Mittelwert & SD für relevante Variablen
for i = 1:length(big5Vars)
    varName = big5Vars{i};
    colIdx = find(strcmp(varNames, varName));
    if ~isempty(colIdx)
        meanRow(colIdx) = mean(summaryTable.(varName), 'omitnan');
        stdRow(colIdx)  = std(summaryTable.(varName), 'omitnan');
    end
end

% --- Neue Tabelle aus einer Zeile erzeugen
meanTable = array2table(meanRow, 'VariableNames', varNames);
stdTable  = array2table(stdRow,  'VariableNames', varNames);

% --- Platzhalter für Teilnehmerinfos
meanTable.VP_Code = "Mean";
stdTable.VP_Code  = "STD";

% --- Zeilen anhängen
summaryTable = [summaryTable; meanTable; stdTable];


% Neue mat-Datei schreiben
save(fullfile(folder, 'questionaireData_all.mat'), 'summaryTable');
disp('Die zusammengefasste mat-Datei wurde erfolgreich erstellt.');
% 
% %% Ausgabeordner für Plots
% plotFolder = fullfile(folder, 'plots');
% if ~exist(plotFolder, 'dir')
%     mkdir(plotFolder);
% end
% 
% % Struktur mit allen Skalen
% factors = struct( ...
%     'Neuroticism', Neuro_Mean, ...
%     'Extraversion', Extra_Mean, ...
%     'Openness', Open_Mean, ...
%     'Agreeableness', Agree_Mean, ...
%     'Conscientiousness', Consc_Mean, ...
%     'OCI', OCI_Sum, ...
%     'Autism', Autism_Sum, ...
%     'GDMS', GDMS_Mean);
% 
% % Alle Histogramme plotten und speichern
% factorNames = fieldnames(factors);
% for i = 1:length(factorNames)
%     name = factorNames{i};
%     values = factors.(name);
% 
%     figure('Visible', 'off');
%     histogram(values, 'BinWidth', 0.5, 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'black');
%     xlabel([name, ' (Mittelwert)']);
%     ylabel('Anzahl der Probanden');
%     title(['Verteilung von ', name]);
%     grid on;
% 
%     % Grafik speichern
%     plotFile = fullfile(plotFolder, [name, '_Histogramm.png']);
%     saveas(gcf, plotFile);
%     close(gcf);
% end
% 
% disp(['Alle Histogramme wurden gespeichert im Ordner: ', plotFolder]);
% 
% %% jittered Boxplots
% % === Boxplot der Big Five mit Einzelpunkten & Fehlerbalken ===
% big5Labels = {'Neuroticism', 'Extraversion', 'Openness', 'Agreeableness', 'Conscientiousness'};
% big5Vars = {'Neuro_Mean', 'Extra_Mean', 'Open_Mean', 'Agree_Mean', 'Consc_Mean'};
% 
% % === Datenmatrix für Boxplots (alle Probanden außer MW & SD Zeilen)
% big5Data = zeros(height(summaryTable)-2, length(big5Vars));
% for i = 1:length(big5Vars)
%     big5Data(:, i) = summaryTable{1:end-2, big5Vars{i}};
% end
% 
% % === Mittelwert & Standardabweichung
% meanVals = mean(big5Data, 'omitnan');
% stdVals = std(big5Data, 'omitnan');
% 
% % === Plot
% figure;
% boxplot(big5Data, 'Labels', big5Labels, 'Colors', 'k', 'Symbol', ''); hold on;
% 
% % Einzelwerte
% n = size(big5Data, 1);
% for i = 1:length(big5Vars)
%     x = repmat(i, n, 1) + (rand(n,1)-0.5)*0.15;
%     scatter(x, big5Data(:, i), 30, 'filled', ...
%         'MarkerFaceColor', [0.2 0.5 0.8], ...
%         'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.6);
% end
% 
% % Fehlerbalken
% errorbar(1:length(big5Vars), meanVals, stdVals, ...
%     'ko', 'MarkerSize', 8, 'LineWidth', 1.5, 'CapSize', 10, 'LineStyle', 'none');
% 
% ylabel('Skalenwert (0–4)');
% title('Big Five: Verteilung, Mittelwert & SD');
% grid on;
% set(gca, 'FontSize', 12);
% 
% % Speichern
% plotFile = fullfile(plotFolder, 'BigFive_Boxplot.png');
% exportgraphics(gcf, plotFile, 'Resolution', 300);
% close(gcf);
% 
% disp(['Big Five Boxplot gespeichert unter: ', plotFile]);
% 
% %% === Normalverteilungsplots für alle Skalen in summaryTable ===
% % Ignoriere Nicht-Skalenvariablen wie VP_Code, Gender, etc.
% excludeVars = {'VP_Code', 'Age', 'Gender', 'Education', 'Hearing', ...
%                'OCI_Binary', 'OCI_Class', 'AutismBinary', 'AutismClass'};
% 
% allVars = summaryTable.Properties.VariableNames;
% numericVars = setdiff(allVars, excludeVars);
% 
% % Nur "echte" Probanden (ohne Mittelwert- und SD-Zeilen)
% isReal = ~ismember(summaryTable.VP_Code, ["Mittelwert", "Standardabweichung"]);
% 
% % Ausgabeordner
% plotFolder = fullfile(folder, 'plots', 'normalverteilung');
% if ~exist(plotFolder, 'dir')
%     mkdir(plotFolder);
% end
% 
% for i = 1:length(numericVars)
%     var = numericVars{i};
%     values = summaryTable{isReal, var};
% 
%     if ~isnumeric(values)
%         continue;
%     end
% 
%     figure('Visible', 'off', 'Name', ['Normalverteilung: ', var]);
% 
%     % 1. Histogramm mit Normalverteilung
%     subplot(1,2,1);
%     histogram(values, 'Normalization', 'pdf', 'FaceColor', [0.6 0.8 1.0]);
%     hold on;
%     x_vals = linspace(min(values), max(values), 100);
%     mu = mean(values, 'omitnan');
%     sigma = std(values, 'omitnan');
%     y_vals = normpdf(x_vals, mu, sigma);
%     plot(x_vals, y_vals, 'r-', 'LineWidth', 2);
%     xlabel(var); ylabel('Dichte');
%     title(['Histogramm mit N(\mu=', num2str(mu, '%.2f'), ', \sigma=', num2str(sigma, '%.2f'), ')']);
% 
%     % 2. Q-Q-Plot
%     subplot(1,2,2);
%     qqplot(values);
%     title('Q-Q-Plot');
%     grid on;
% 
%     % Speichern
%     saveas(gcf, fullfile(plotFolder, [var '_Normalverteilung.png']));
%     close(gcf);
% end
% 
% disp(['Alle Normalverteilungsplots wurden gespeichert unter: ', plotFolder]);
