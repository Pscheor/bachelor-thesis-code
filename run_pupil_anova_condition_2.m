function run_pupil_anova_condition_2(trial)
% Mixed rm ANOVA (no binning) on Pupil_base_avg
% Within:  Condition (Liberal, Conservative)
% Between: ExtraversionGroup2 (introvert, extravert)  [ambiverts excluded]
% Plots: Liberal boxplot | Conservative boxplot | Interaction (side-by-side)

    fprintf('\n==== Mixed ANOVA (no binning): Pupil size by Condition x ExtraversionGroup ====\n');

    pupilVar = 'Pupil_base_avg';
    if ~ismember(pupilVar, trial.Properties.VariableNames)
        error('Variable "%s" not found in table.', pupilVar);
    end

    %% Filter (no ambiverts) & normalize
    grp  = lower(strtrim(string(trial.ExtraversionGroup2)));
    cond = lower(strtrim(string(trial.Condition)));
    keep = ismember(grp, {'introvert','extravert'}) & ismember(cond, {'liberal','conservative'});
    trial = trial(keep, :);
    if isempty(trial)
        error('No data after filtering (introvert/extravert & liberal/conservative).');
    end

    trial.ExtraversionGroup2 = categorical(lower(strtrim(string(trial.ExtraversionGroup2))));
    trial.Condition          = categorical(lower(strtrim(string(trial.Condition))));
    trial.Condition          = renamecats(trial.Condition, {'liberal','conservative'}, {'Liberal','Conservative'});

    %% Aggregate: mean pupil per VP & condition
    VPs = unique(trial.VPCode);
    wide = table();
    for i = 1:numel(VPs)
        vp  = VPs(i);
        vpD = trial(trial.VPCode == vp, :);

        g = unique(vpD.ExtraversionGroup2);
        if numel(g) ~= 1 || ismissing(g)
            continue;
        end

        vL = vpD.(pupilVar)(vpD.Condition=='Liberal');
        vC = vpD.(pupilVar)(vpD.Condition=='Conservative');
        if isempty(vL) || isempty(vC)
            continue;
        end

        mL = mean(vL,'omitnan');
        mC = mean(vC,'omitnan');
        if ~isfinite(mL) || ~isfinite(mC)
            continue;
        end

        wide = [wide; table(string(vp), categorical(g), mL, mC, ...
            'VariableNames', {'VPCode','Group','Liberal','Conservative'})]; %#ok<AGROW>
    end

    if isempty(wide)
        error('No valid participants for the Mixed ANOVA.');
    end

    % order groups
    desiredOrder = intersect({'introvert','extravert'}, cellstr(categories(wide.Group)), 'stable');
    wide.Group = reordercats(wide.Group, desiredOrder);

    %% Mixed ANOVA
    within = table(categorical({'Liberal';'Conservative'}), 'VariableNames', {'Condition'});
    rm = fitrm(wide, 'Liberal,Conservative ~ Group', 'WithinDesign', within);
    ranovatbl = ranova(rm, 'WithinModel', 'Condition');

    fprintf('\n--- Mixed ANOVA table ---\n');
    disp(ranovatbl);

    %% Partial eta^2
    fprintf('\n--- Effect sizes (partial eta^2) ---\n');
    rnames = ranovatbl.Properties.RowNames;
    for i = 1:height(ranovatbl)
        row = rnames{i};
        if contains(row,'Group') || strcmp(row,'Condition') || contains(row,':')
            SS_eff = ranovatbl.SumSq(i);
            if contains(row,':') || strcmp(row,'Condition')
                SS_err = ranovatbl.SumSq(strcmp(rnames,'Error(Condition)'));
            else
                SS_err = ranovatbl.SumSq(strcmp(rnames,'Error'));
            end
            eta2 = SS_eff / (SS_eff + SS_err);
            fprintf('%-18s partial eta^2 = %.4f\n', row, eta2);
        end
    end

    %% Post-hoc tests
    fprintf('\n--- Post-hoc: paired t-tests within groups (Liberal vs Conservative) ---\n');
    grCats = categories(wide.Group);
    for g = 1:numel(grCats)
        gname = grCats{g};
        idx = (wide.Group == gname);
        [~, p, ~, st] = ttest(wide.Liberal(idx), wide.Conservative(idx));
        diffs = wide.Liberal(idx) - wide.Conservative(idx);
        dz = mean(diffs,'omitnan') / std(diffs,'omitnan');
        fprintf('%s: t(%d) = %.3f, p = %.4f, Cohen''s dz = %.3f\n', gname, st.df, st.tstat, p, dz);
    end

    fprintf('\n--- Post-hoc: between-group comparisons per condition (Welch t-test) ---\n');
    if numel(grCats) >= 2
        g1 = grCats{1}; g2 = grCats{2};
        L1 = wide.Liberal(wide.Group==g1); L2 = wide.Liberal(wide.Group==g2);
        [~, pL, ~, stL] = ttest2(L1, L2, 'Vartype','unequal');
        dL = cohens_d(L1, L2);
        fprintf('Liberal:     t(%d) = %.3f, p = %.4f, Cohen''s d = %.3f\n', round(stL.df), stL.tstat, pL, dL);

        C1 = wide.Conservative(wide.Group==g1); C2 = wide.Conservative(wide.Group==g2);
        [~, pC, ~, stC] = ttest2(C1, C2, 'Vartype','unequal');
        dC = cohens_d(C1, C2);
        fprintf('Conservative: t(%d) = %.3f, p = %.4f, Cohen''s d = %.3f\n', round(stC.df), stC.tstat, pC, dC);
    end

    %% ---- PLOTS (1 x 3 side-by-side) ----
    colorIntro = [0.15 0.55 0.85];
    colorExtra = [0.90 0.40 0.25];
    colors = [colorIntro; colorExtra];

    % Hier definierst du die Labels manuell
    groupLabels = {'Introvert','Extravert'};

    figure('Name','Pupil size by Condition x Extraversion', ...
           'Units','normalized','Position',[0.05 0.35 0.90 0.45]);
    tl = tiledlayout(1,3, 'TileSpacing','compact','Padding','compact');
    title(tl, 'Pupil size by Condition x Extraversion Group');

    % Liberal
    ax1 = nexttile(tl,1);
    valsL_intro = wide.Liberal(wide.Group=='introvert');
    valsL_extra = wide.Liberal(wide.Group=='extravert');
    plot_box_two_groups(ax1, valsL_intro, valsL_extra, colors, 'Liberal', 'Pupil size (mean)', groupLabels);

    % Conservative
    ax2 = nexttile(tl,2);
    valsC_intro = wide.Conservative(wide.Group=='introvert');
    valsC_extra = wide.Conservative(wide.Group=='extravert');
    plot_box_two_groups(ax2, valsC_intro, valsC_extra, colors, 'Conservative', 'Pupil size (mean)', groupLabels);

    % Interaction
    ax3 = nexttile(tl,3); 
    hold(ax3,'on');
    mIntro = [mean(valsL_intro,'omitnan'), mean(valsC_intro,'omitnan')];
    mExtra = [mean(valsL_extra,'omitnan'), mean(valsC_extra,'omitnan')];
    sIntro = [std(valsL_intro,'omitnan')/sqrt(max(1,sum(~isnan(valsL_intro)))), ...
              std(valsC_intro,'omitnan')/sqrt(max(1,sum(~isnan(valsC_intro))))];
    sExtra = [std(valsL_extra,'omitnan')/sqrt(max(1,sum(~isnan(valsL_extra)))), ...
              std(valsC_extra,'omitnan')/sqrt(max(1,sum(~isnan(valsC_extra))))];
    xpts = [1 2];
    errorbar(ax3, xpts, mIntro, sIntro, '-o', 'Color', colorIntro, ...
             'MarkerFaceColor', colorIntro, 'LineWidth', 1.8);
    errorbar(ax3, xpts, mExtra, sExtra, '-o', 'Color', colorExtra, ...
             'MarkerFaceColor', colorExtra, 'LineWidth', 1.8);
    set(ax3,'XTick',xpts,'XTickLabel',{'Liberal','Conservative'});
    ylabel(ax3,'Pupil size (mean)');
    title(ax3,'Interaction (Mean ± SEM)');
    grid(ax3,'on'); 
    box(ax3,'on');
    legend(ax3, groupLabels, 'Location','best');

    % Set specific axis limits
    axis(ax3, [0.5, 2.5, 4550, 5000]);

end

%% ===== Lokale Hilfsfunktionen =====
function d = cohens_d(a, b)
    a = a(:); b = b(:);
    na = sum(~isnan(a)); nb = sum(~isnan(b));
    sa2 = var(a,'omitnan'); sb2 = var(b,'omitnan');
    sp = sqrt(((na-1)*sa2 + (nb-1)*sb2) / max(1,(na+nb-2)));
    d = (mean(a,'omitnan') - mean(b,'omitnan')) / sp;
end

function plot_box_two_groups(ax, valsIntro, valsExtra, colors, condName, ylab, xLabels)
    axes(ax); cla(ax); hold(ax,'on');
    y  = [valsIntro(:); valsExtra(:)];
    gx = [repmat({xLabels{1}}, numel(valsIntro),1); ...
          repmat({xLabels{2}}, numel(valsExtra),1)];
    boxplot(y, gx, 'GroupOrder', xLabels, ...
            'Colors',[0 0 0], 'Symbol','k+', 'Widths',0.6);
    h = findobj(ax,'Tag','Box');
    for j = 1:min(numel(h),2)
        patch(get(h(j),'XData'), get(h(j),'YData'), colors(3-j,:), ...
              'FaceAlpha',0.25, 'EdgeColor', colors(3-j,:));
    end
    x = [ones(numel(valsIntro),1); 2*ones(numel(valsExtra),1)];
    jitter = (rand(size(x))-0.5)*0.25;
    Cpts = colors(x,:);
    scatter(x+jitter, y, 36, Cpts, 'filled', 'MarkerFaceAlpha',0.65, ...
            'MarkerEdgeColor',[0 0 0], 'MarkerEdgeAlpha',0.15);

    yl = ylim(ax); yr = diff(yl); if yr==0, yr=1; end
    for k = 1:2
        if k==1, vals = valsIntro; else, vals = valsExtra; end
        mu = mean(vals,'omitnan');
        sem = std(vals,'omitnan')/max(1,sqrt(sum(~isnan(vals))));
        plot(k, mu, 'kd', 'MarkerFaceColor', colors(k,:), 'MarkerSize', 8);
        errorbar(k, mu, sem, 'k', 'LineWidth', 1.1, 'CapSize', 10);
        text(k, yl(1)+0.05*yr, sprintf('n=%d', sum(~isnan(vals))), ...
            'HorizontalAlignment','center','Color',[0.25 0.25 0.25], 'FontSize',9);
    end

    set(ax,'XTick',[1 2],'XTickLabel',xLabels);
    ylabel(ax, ylab);
    title(ax, condName);
    grid(ax,'on'); box(ax,'on');
    ylim(ax, [yl(1)-0.06*yr, yl(2)+0.08*yr]);
    hold(ax,'off');
end
