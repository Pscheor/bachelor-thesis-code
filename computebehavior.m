function [behavior_block, behavior_condition, behavior_trial] = computebehavior(cfg)

% === Daten laden
load(cfg.datafile);  % erwartet: Data, Param
VPCode = string(Param.SubjInfo.SubjectCode);
PPNumber = extractBetween(cfg.datafile, 'PP_', filesep);  % z.B. "001"

% === Blockfilterung: nur Blöcke mit >= 50 Trials behalten
for t = 1:numel(Data.Task)
    res = Data.Task(t).Results;
    flow = isfield(Data.Task(t), 'Flow') && ~isempty(Data.Task(t).Flow);
    flowData = [];
    if flow
        flowData = Data.Task(t).Flow;
    end

    isValid = cellfun(@(tab) size(tab,1) >= 50, res);
    res = res(isValid);
    flowData = flowData(isValid);

    switch t
        case 1
            if numel(res) > 1
                res = res(end);
                flowData = flowData(end);
            end
        case 2
            if numel(res) > 4
                res = res(end-3:end);
                flowData = flowData(end-3:end);
            end
        case 3
            if numel(res) > 4
                res = res(end-3:end);
                flowData = flowData(end-3:end);
            end
    end

    Data.Task(t).Results = res;
    Data.Task(t).Flow = flowData;
end

% === SDT-Codes berechnen & Bereinigung
for t = 1:numel(Data.Task)
    for r = 1:numel(Data.Task(t).Results)
        tab = Data.Task(t).Results{r};

        % --- Erste 2 Trials entfernen
        if size(tab,1) > 2
            tab = tab(3:end,:);
        else
            tab = [];  % zu wenig Daten
        end

        if isempty(tab)
            Data.Task(t).Results{r} = [];
            continue
        end

        stim = tab(:,2);     % 1 = Ton anwesend, 0 = nicht
        resp = tab(:,3);     % 1 = Ja, 2 = Nein, 0 = keine Antwort

        code = nan(size(stim));
        code(resp == 0) = 4;
        code(resp == 1 & stim == 1) = 0;
        code(resp == 1 & stim == 0) = 2;
        code(resp == 2 & stim == 1) = 1;
        code(resp == 2 & stim == 0) = 3;

        tab(:,10) = code;

        tab(tab(:,6) < 0, 6) = NaN;
        tab(resp == 0, 4)    = NaN;

        Data.Task(t).Results{r} = tab;
    end
end

% === Einzelblöcke: behavior_block erstellen
behavior_block = table();
for t = 1:numel(Data.Task)
    if isempty(Data.Task(t).Results), continue; end
    switch t
        case 2, cond = "Liberal";
        case 3, cond = "Conservative";
        otherwise, continue;
    end

    for r = 1:numel(Data.Task(t).Results)
        tab = Data.Task(t).Results{r};
        if isempty(tab), continue; end
        rt = tab(:,4); rt(rt < 0.01) = NaN;
        [dprime, c, nHit, nMiss, nCR, nFA] = computeSDT(tab(:,10));

        row = table();
        row.Filename = "T" + t + "_R" + r;
        row.Condition = cond;
        row.VPCode = VPCode;
        row.MeanRT = nanmean(rt);
        row.DimmFactor = nanmean(tab(:,5));
        row.MeanConfidence = nanmean(tab(:,6));
        row.MeanReactionTime = nanmean(tab(:,7));
        row.dPrime = dprime;
        row.c = c;

        if r <= numel(Data.Task(t).Flow) && isfield(Data.Task(t).Flow(r), 'Rating')
            rating = double(Data.Task(t).Flow(r).Rating);
            if isnumeric(rating) && rating >= 0 && rating <= 100
                row.Flow = rating;
            else
                row.Flow = NaN;
            end
        else
            row.Flow = NaN;
        end

        row.nHit = nHit; row.nMiss = nMiss;
        row.nCR = nCR; row.nFA = nFA;

        behavior_block = [behavior_block; row];
    end
end

% === Condition Table: eine Zeile pro VP
condLabels = ["Liberal", "Conservative"];
conditionRow = table();
conditionRow.VPCode = VPCode;

for i = 1:2
    cond = condLabels(i);
    sub = behavior_block(strcmp(behavior_block.Condition, cond), :);
    if isempty(sub), continue; end

    suffix = lower(cond);
    conditionRow.("MeanRT_" + suffix) = nanmean(sub.MeanRT);
    conditionRow.("DimmFactor_" + suffix) = nanmean(sub.DimmFactor);
    conditionRow.("MeanConfidence_" + suffix) = nanmean(sub.MeanConfidence);
    conditionRow.("MeanReactionTime_" + suffix) = nanmean(sub.MeanReactionTime);
    conditionRow.("dPrime_" + suffix) = nanmean(sub.dPrime);
    conditionRow.("c_" + suffix) = nanmean(sub.c);
    conditionRow.("Flow_" + suffix) = nanmean(sub.Flow);
end

behavior_condition = conditionRow;

% === Trial-Tabelle: behavior_trial erstellen
behavior_trial = table();

for t = 1:numel(Data.Task)
    if ~ismember(t, [2, 3]), continue; end  % Nur Liberal & Conservative

    condition = ["Liberal", "Conservative"];
    condLabel = condition(t-1);

    res = Data.Task(t).Results;
    flowData = Data.Task(t).Flow;

    for r = 1:numel(res)
        blockTab = res{r};
        blockNum = r;

        stim = blockTab(:,2);       % Stimulus
        resp = blockTab(:,3);       % Response
        rt   = blockTab(:,4);       % RT
        conf = blockTab(:,6);       % Confidence

        conf(conf < 0) = NaN;
        rt(resp == 0) = NaN;

        rating = NaN;
        if isfield(flowData(r), 'Rating')
            tempRating = flowData(r).Rating;
            if isnumeric(tempRating) && tempRating >= 0 && tempRating <= 100
                rating = tempRating;
            end
        end

        for tr = 1:size(stim,1)
            newRow = table();
            newRow.VPCode = VPCode;
            newRow.PP_ID = string(PPNumber);
            newRow.Condition = condLabel;
            newRow.BlockNumber = blockNum;
            newRow.TrialNumber = tr;
            newRow.Stimulus = stim(tr);
            newRow.Response = resp(tr);
            newRow.ReactionTime = rt(tr);
            newRow.Confidence = conf(tr);
            newRow.DimmFactor = blockTab(tr,5);
            newRow.Flow = rating;
        
            % Treffer-Korrektheit
            if resp(tr) == 1 && stim(tr) == 1
                newRow.Correct = 1;
            elseif resp(tr) == 2 && stim(tr) == 0
                newRow.Correct = 1;
            elseif resp(tr) == 0
                newRow.Correct = NaN;
            else
                newRow.Correct = 0;
            end
        
            % === SDT-Code berechnen
            code = NaN;
            if resp(tr) == 0
                code = 4;
            elseif resp(tr) == 1 && stim(tr) == 1
                code = 0;
            elseif resp(tr) == 2 && stim(tr) == 1
                code = 1;
            elseif resp(tr) == 1 && stim(tr) == 0
                code = 2;
            elseif resp(tr) == 2 && stim(tr) == 0
                code = 3;
            end
            newRow.SDTCode = code;
        
            behavior_trial = [behavior_trial; newRow];
        end

    end
end

% === Vorherige Daten laden (optional)
fileBlockData = fullfile(cfg.datapath, 'behavior_block.mat');
fileConditionData = fullfile(cfg.datapath, 'behavior_condition.mat');
fileTrialData = fullfile(cfg.datapath, 'behavior_trial.mat');

if isfile(fileBlockData)
    S = load(fileBlockData, 'behavior_block');
    behavior_block = [S.behavior_block; behavior_block];
end

if isfile(fileConditionData)
    S = load(fileConditionData, 'behavior_condition');
    behavior_condition = [S.behavior_condition; behavior_condition];
end

if isfile(fileTrialData)
    S = load(fileTrialData, 'behavior_trial');
    behavior_trial = [S.behavior_trial; behavior_trial];
end

% === Speichern
try
    save(fileBlockData, 'behavior_block');
    save(fileConditionData, 'behavior_condition');
    save(fileTrialData, 'behavior_trial');
    fprintf("✅ Tabellen gespeichert: behavior_block, behavior_condition, behavior_trial\n");
catch ME
    warning("❌ Fehler beim Speichern: %s", ME.message);
end

end
