clear 
restoredefaultpath
% Set default interpreter for title, xlabel, ylabel, zlabel
set(groot, 'defaultAxesTickLabelInterpreter', 'none');  % For tick labels
set(groot, 'defaultTextInterpreter', 'none');  % For text and title
set(groot, 'defaultLegendInterpreter', 'none');  % For legend

if ispc  
  datapath = "C:\Uni\BA\Daten";
  toolspath = "C:\MATLAB\tools";
  % addpath(fullfile(toolspath, 'fieldtrip-20231025')); ft_defaults
  addpath(fullfile(toolspath, 'fieldtrip-20240603')); ft_defaults
  plotpath = "C:\Uni\BA\Plots";
elseif ismac
  datapath = "/Users/kloosterman/projectdata/pupilperso";
  toolspath = "/Users/kloosterman/Documents/GitHub";
  % addpath(fullfile(toolspath, 'fieldtrip-20231025')); ft_defaults
  addpath(fullfile(toolspath, 'FieldTrip')); ft_defaults  
  plotpath = '/Users/kloosterman/Library/CloudStorage/Dropbox/PROJECTS';
end

addpath(genpath(fullfile(toolspath, 'PupilPerso')));

% define subjects: enter number here
SUBJ = {};
for i = 31
    subjID = sprintf('%03d', i);
    SUBJ{end+1} = subjID;  % ✅ wichtig!
end

SUBJbool = true(size(SUBJ));
SUBJ = SUBJ(SUBJbool);

disp(SUBJ)

%% preprocessing, compute entropy, freqanalysis and ERPs, or computebehavior
% Which processing step?
fun2run = @preprocessing; %  computemixedmodel preprocessing computemMSE computebehavior merge_subjects
overwrite = 1;
cfg = [];
cfg.datapath = datapath;
cfglist = {};
for isub = 1:length(SUBJ)
  cfg.subjno = SUBJ{isub};
  cfg.outfile = sprintf('SUB%s_cond', cfg.subjno);
  % cfg.datafile = fullfile(datapath, 'data', 'pilot1', SUBJ{isub}, sprintf('PupilPersonality_SubjID%s', SUBJ{isub}));
  % cfg.datafile = fullfile(datapath, 'data', 'pilot2', SUBJ{isub}, sprintf('PupilPersonality_SubjID%s', SUBJ{isub}));
  % === Robuste Suche nach MAT-Dateien (mit oder ohne führende Nullen)
  subjectFolder = fullfile(datapath, ['PP_', SUBJ{isub}]); % input
  outputFolder = fullfile(datapath, 'preproc'); % output
  mkdir(outputFolder);

  % Suche nach allen Dateien, die SubjID enthalten – egal wie viele Nullen
  filematch = dir(fullfile(subjectFolder, sprintf('PupilPersonality_SubjID*%s*.mat', SUBJ{isub})));

  % Wenn keine gefunden wurde, alternativ nochmal mit exakt SUBJ{isub}
  if isempty(filematch)
    filematch = dir(fullfile(subjectFolder, sprintf('PupilPersonality_SubjID%s_*.mat', SUBJ{isub})));
  end

  if isempty(filematch)
    warning('Keine passende .mat-Datei für VP %s gefunden.', SUBJ{isub});
    continue;
  end

  cfg.datafile = fullfile(filematch(1).folder, filematch(1).name);
  cfg.outfile = sprintf('behavior_%s', cfg.subjno);
  cfg.outputFolder = outputFolder;

  if overwrite;
    cfglist{end+1} = cfg
  end % TODO if ismerge; load data; else cfglist{end+1} = cfg
end
% cfglist = cfglist(1:3)

if strcmp(func2str(fun2run), 'merge_subjects')
  eeg = fun2run(cfglist);
  behav = eeg.mse.across.AA.trialinfo;
elseif strcmp(func2str(fun2run), 'computebehavior')
  % Ergebnisse und Blockinfos separat sammeln
  behav = cell(size(cfglist));
  blockinfo_all = cell(size(cfglist));

  for i = 1:numel(cfglist)
    [behav{i}, blockinfo_all{i}] = fun2run(cfglist{i});
  end

  behav = vertcat(behav{:});
  blockinfo_all = vertcat(blockinfo_all{:});

  % Gemeinsame Blockinfo speichern
  save(fullfile(datapath, 'blockinfo.mat'), 'blockinfo_all');
  fprintf('📄 Gesamte Blockinfo gespeichert: blockinfo.mat\n');

else
  trialresults = cellfun(fun2run, cfglist, 'UniformOutput', false);
  trialresults = vertcat(trialresults{:});
end

