function results = analyze_betas_group_condition(trial)
% Analyse von linearen & quadratischen Betas (aus 6 Bins) pro VP × Condition (Lib/Con)
% Verwendet nur introvert/extravert aus ExtraversionGroup2.
% Gibt Tabelle 'results' zurück und erzeugt Boxplots & Linienplots.
%
% Voraussetzungen:
% - trial: Tabelle mit VPCode, ExtraversionGroup2, Condition, Pupil_base_avg, SDTCode
% - computeSDT(codes) vorhanden

fprintf('\n==== Starte Analyse der Beta-Koeffizienten (Group2 × Condition) ====\n');

%% -------------------- Filter & Vorbereitung --------------------
% nur introvert/extravert behalten
maskValid = ismember(lower(strtrim(string(trial.ExtraversionGroup2))), ["introvert","extravert"]);
trial = trial(maskValid, :);

% Vereinheitlichen
trial.ExtraversionGroup2 = categorical(lower(strtrim(string(trial.ExtraversionGroup2))));
trial.Condition          = categorical(strtrim(string(trial.Condition)));
trial = trial(ismember(trial.Condition, categorical({'Liberal','Conservative'})), :);

numBins = 6;
latVar  = 'Pupil_base_avg';
conds   = {'Liberal','Conservative'};
groups  = {'introvert','extravert'};

% Ergebnis-Tabelle robust initialisieren (damit Property-Zugriffe sicher sind)
results = table('Size',[0 6], ...
    'VariableTypes', {'string','categorical','categorical','double','double','double'}, ...
    'VariableNames', {'VPCode','Group','Condition','Beta_Linear','Beta_Quad1','Beta_Quad2'});

% Farben
col.introvert_Liberal      = [0.6 0.8 1.0];  % hellblau
col.introvert_Conservative = [0.0 0.2 0.6];  % dunkelblau
col.extravert_Liberal      = [1.0 0.6 0.2];  % orange
col.extravert_Conservative = [0.8 0.1 0.1];  % rot

%% -------------------- Betas pro VP × Condition --------------------
VPs = unique(trial.VPCode);
for i = 1:numel(VPs)
    vp = VPs(i);
    vpData = trial(strcmp(trial.VPCode, vp), :);

    % Gruppen-Kategorie pro VP (ist bereits categorical in der Tabelle)
    gvals = unique(vpData.ExtraversionGroup2);
    gvals = gvals(~ismember(gvals, categorical(missing)));
    if numel(gvals) ~= 1
        % uneindeutig / fehlend -> überspringen
        continue
    end
    gCat = gvals(1);  % categorical

    for c = 1:numel(conds)
        condName = conds{c};
        condMask = (vpData.Condition == condName);
        condData = vpData(condMask, :);

        if height(condData) < numBins
            continue
        end

        pupil = condData.(latVar);
        codes = condData.SDTCode;

        if all(isnan(pupil)) || numel(pupil) < numBins
            continue
        end

        % Binning auf Quantile (robust mit Fallback)
        try
            edges = [-inf, quantile(pupil, (1:numBins-1)/numBins), inf];
        catch
            q = prctile(pupil, (1:(numBins-1))/numBins*100);
            edges = [-inf, q, inf];
        end
        binIdx = discretize(pupil, edges);

        pupilMeans = NaN(1,numBins);
        dprimes    = NaN(1,numBins);
        for b = 1:numBins
            idx = (binIdx == b);
            if any(idx)
                pupilMeans(b) = mean(pupil(idx), 'omitnan');
                [dprimes(b), ~] = computeSDT(codes(idx));
            end
        end

        if any(isnan(pupilMeans)) || any(isnan(dprimes))
            continue
        end

        % Zentrieren
        pupil_dm = pupilMeans - mean(pupilMeans, 'omitnan');

        % Lineares Modell
        mdl_lin  = fitlm(pupil_dm, dprimes, 'linear');
        beta_lin = mdl_lin.Coefficients.Estimate(2);

        % Quadratisches Modell
        mdl_quad = fitlm(pupil_dm, dprimes, 'quadratic');
        beta_q1  = mdl_quad.Coefficients.Estimate(2);
        beta_q2  = mdl_quad.Coefficients.Estimate(3);

        % neue Zeile als TABLE (keine Cell) mit korrekten Typen
        newRow = table(string(vp), gCat, categorical(string(condName)), ...
                       beta_lin, beta_q1, beta_q2, ...
                       'VariableNames', results.Properties.VariableNames);

        results = [results; newRow]; %#ok<AGROW>
    end
end

if isempty(results)
    error('Keine gültigen Betas berechnet. Prüfe Spalten (Condition, Pupil_base_avg, SDTCode) und Datenumfang.');
end
fprintf('✅ %d Zeilen in results (VP × Condition) erstellt.\n', height(results));

%% -------------------- Statistik: T-Tests & Effekte --------------------
betaTypes = {'Beta_Linear','Beta_Quad1','Beta_Quad2'};
fprintf('\n==== Statistische Tests ====\n');

for b = 1:numel(betaTypes)
    bname = betaTypes{b};
    fprintf('\n--- %s ---\n', bname);

    % Gruppenvergleich intro vs extravert je Condition (Welch)
    for c = 1:numel(conds)
        condName = conds{c};
        vals_intro = results.(bname)(results.Group=='introvert' & results.Condition==condName);
        vals_extra = results.(bname)(results.Group=='extravert' & results.Condition==condName);

        if sum(~isnan(vals_intro)) >= 2 && sum(~isnan(vals_extra)) >= 2
            [~, p, ~, st] = ttest2(vals_intro, vals_extra, 'Vartype','unequal');
            d = cohens_d_two(vals_intro, vals_extra);
            fprintf('Intro vs Extra (%s): t(%.1f) = %.3f, p = %.4f, d = %.3f\n', ...
                condName, st.df, st.tstat, p, d);
        else
            fprintf('Intro vs Extra (%s): zu wenige Daten.\n', condName);
        end
    end

    % One-sample Tests gegen 0 je (Group × Condition)
    for gi = 1:numel(groups)
        gname = groups{gi};
        for c = 1:numel(conds)
            condName = conds{c};
            data = results.(bname)(results.Group==gname & results.Condition==condName);
            n = sum(~isnan(data));
            if n >= 2
                [~, p0, ~, st0] = ttest(data);
                eta2 = (st0.tstat^2) / (st0.tstat^2 + st0.df);
                d0   = mean(data,'omitnan') / std(data,'omitnan'); % Cohen's d (one-sample)
                fprintf('Gegen 0: %s – %s: t(%d) = %.3f, p = %.4f, d = %.3f, η² = %.3f\n', ...
                    gname, condName, st0.df, st0.tstat, p0, d0, eta2);
            else
                fprintf('Gegen 0: %s – %s: zu wenige Daten.\n', gname, condName);
            end
        end
    end
end

%% -------------------- Boxplots --------------------
figure('Name','Boxplots: Beta-Koeffizienten','Units','normalized','Position',[0.08 0.42 0.84 0.48]);
for b = 1:numel(betaTypes)
    subplot(1,3,b); hold on; grid on; box on
    bname = betaTypes{b};
    title(strrep(bname,'_','\_')); ylabel('Beta');

    defs = {
        'introvert','Liberal',      col.introvert_Liberal,      1;
        'introvert','Conservative', col.introvert_Conservative, 2;
        'extravert','Liberal',      col.extravert_Liberal,      3;
        'extravert','Conservative', col.extravert_Conservative, 4;
    };

    for k = 1:size(defs,1)
        gname = defs{k,1};
        cname = defs{k,2};
        clr   = defs{k,3};
        pos   = defs{k,4};

        vals = results.(bname)(results.Group==gname & results.Condition==cname);
        vals = vals(~isnan(vals));
        if isempty(vals), continue; end

        % Boxchart: ein Boxplot an Position pos
        x = pos * ones(size(vals));
        boxchart(x, vals, 'BoxFaceColor', clr, 'WhiskerLineColor','k', 'BoxEdgeColor','k');
        % leichte Punktwolke
        jitter = (rand(size(vals))-0.5)*0.12;
        scatter(x + jitter, vals, 24, clr, 'filled', 'MarkerFaceAlpha',0.65, 'MarkerEdgeColor','k', 'MarkerEdgeAlpha',0.15);
    end

    xlim([0.5 4.5]);
    set(gca,'XTick',1:4,'XTickLabel',{'Intro-Lib','Intro-Con','Extra-Lib','Extra-Con'},'XTickLabelRotation',35);
end

%% -------------------- Linienplots (Mittel ± SEM) --------------------
figure('Name','Beta Mittelwerte ± SEM','Units','normalized','Position',[0.12 0.18 0.8 0.42]);
for b = 1:numel(betaTypes)
    subplot(1,3,b); hold on; grid on; box on
    bname = betaTypes{b};
    title(strrep(bname,'_','\_')); ylabel('Beta');

    % Introvert
    vIL = results.(bname)(results.Group=='introvert' & results.Condition=='Liberal');
    vIC = results.(bname)(results.Group=='introvert' & results.Condition=='Conservative');
    mI = [mean(vIL,'omitnan'), mean(vIC,'omitnan')];
    sI = [std(vIL,'omitnan')/max(1,sqrt(sum(~isnan(vIL)))), ...
          std(vIC,'omitnan')/max(1,sqrt(sum(~isnan(vIC))))];
    plot(1:2, mI, '-o', 'Color',[0.2 0.2 0.2],'LineWidth',1.6,'MarkerFaceColor',col.introvert_Liberal,'MarkerEdgeColor','k');
    errorbar(1:2, mI, sI, 'k','LineStyle','none','CapSize',8);
    scatter(2, mI(2), 60, col.introvert_Conservative, 'filled', 'MarkerEdgeColor','k');

    % Extravert
    vEL = results.(bname)(results.Group=='extravert' & results.Condition=='Liberal');
    vEC = results.(bname)(results.Group=='extravert' & results.Condition=='Conservative');
    mE = [mean(vEL,'omitnan'), mean(vEC,'omitnan')];
    sE = [std(vEL,'omitnan')/max(1,sqrt(sum(~isnan(vEL)))), ...
          std(vEC,'omitnan')/max(1,sqrt(sum(~isnan(vEC))))];
    plot(1:2, mE, '-o', 'Color',[0.2 0.2 0.2],'LineWidth',1.6,'MarkerFaceColor',col.extravert_Liberal,'MarkerEdgeColor','k');
    errorbar(1:2, mE, sE, 'k','LineStyle','none','CapSize',8);
    scatter(2, mE(2), 60, col.extravert_Conservative, 'filled', 'MarkerEdgeColor','k');

    set(gca,'XTick',1:2,'XTickLabel',{'Liberal','Conservative'});
    legend({'Introvert','Extravert'}, 'Location','best');
end

fprintf('\n==== Fertig. "results" enthält die Betas; Plots erstellt. ====\n');
end

%% ================= Hilfsfunktion =================
function d = cohens_d_two(a, b)
% Cohen's d (unabhängige Stichproben)
    a = a(:); b = b(:);
    a = a(~isnan(a)); b = b(~isnan(b));
    na = numel(a); nb = numel(b);
    if na<2 || nb<2
        d = NaN; return
    end
    sa2 = var(a, 1); % Pop-Variante für robuste sp-Kombination
    sb2 = var(b, 1);
    sp  = sqrt(((na-1)*sa2 + (nb-1)*sb2) / max(1,(na+nb-2)));
    d   = (mean(a) - mean(b)) / sp;
end
