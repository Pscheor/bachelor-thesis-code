function [data] = preprocessing(cfg)
% load EyeLink edf data into FieldTrip, interpolate blinks and filter

% path to edf2asc tool to convert edf to asc
if ispc
  edf2asc = 'C:\MATLAB\tools\PupilPerso\edf2asc.exe';
elseif ismac
  edf2asc = '/Users/kloosterman/Library/CloudStorage/Dropbox/MATLAB/toolbox/ETanalysis/edfapi_universal/edf2asc';
end

% unpack variables contained in input cfg
outfile = cfg.outfile;
outputFolder = cfg.outputFolder;

datafile = cfg.datafile;
subjno = cfg.subjno;
disp(subjno)

datafolder = fileparts(datafile);  % Speichere Ordnerpfad
cd(datafolder);                    % Wechsle nur einmal in den Ordner

%% SubjCode aus zentraler Datei extrahieren
infofile = dir(fullfile(fileparts(datafile), 'PupilPersonality_SubjID*_????.mat'));

if isempty(infofile)
   error('Keine passende Meta-Datei mit SubjectInfo im Ordner %s gefunden.', fileparts(datafile));
elseif numel(infofile) > 1
   warning('Mehrere passende .mat-Dateien gefunden, erste wird verwendet: %s', infofile(1).name);
end

% Lade die Datei mit Param.SubjInfo.SubjectCode
loadedInfo = load(fullfile(infofile(1).folder, infofile(1).name));

if ~isfield(loadedInfo, 'Param') || ~isfield(loadedInfo.Param, 'SubjInfo') || ~isfield(loadedInfo.Param.SubjInfo, 'SubjectCode')
    error('SubjectCode konnte in %s nicht gefunden werden.', infofile(1).name);
end

SubjectCode = loadedInfo.Param.SubjInfo.SubjectCode;
fprintf('SubjectCode extrahiert: %s\n', SubjectCode);

if isempty(infofile)
end

%% edf2asc
conds = ["liberal" "Konservative" "Baseline" "Training"];
data_avg = {};
timelock_avg = {};

for icond = 1:3
  edflist = dir( "*" + conds(icond) + "*.edf");

  alldata = {};
  all_timelock = {};

  for irun = 1:length(edflist)
    % for irun = 2:3
    [~, eyename] = fileparts(edflist(irun).name);
    filename_asc = [eyename '.asc'];

    % Check if ASC file already exists
    if ~isfile(filename_asc)
      fprintf('Konvertiere %s nach %s...\n', edflist(irun).name, filename_asc);
      % Starte Konvertierung
      system(sprintf('"%s" -y "%s"', edf2asc, edflist(irun).name));

      % Warte darauf, dass die ASC-Datei erstellt wurde
      max_wait_time = 10; % in Sekunden
      elapsed = 0;
      while ~isfile(filename_asc) && elapsed < max_wait_time
        pause(0.5);
        elapsed = elapsed + 0.5;
      end

      % Validierung: Datei da und >0 Byte?
      if ~isfile(filename_asc) || dir(filename_asc).bytes < 1000
        warning('Konvertierung von %s fehlgeschlagen oder Datei unvollständig.', edflist(irun).name);
        continue;
      end

      % Warte bis ASC-Datei fertig erstellt ist
      max_wait = 5; % max 5 Sekunden warten
      wait_time = 0;
      while ~isfile(filename_asc) && wait_time < max_wait
        pause(0.5); % 0.5 Sekunde warten
        wait_time = wait_time + 0.5;
      end

      if ~isfile(filename_asc)
        warning('ASC-Datei %s wurde nicht rechtzeitig erstellt – überspringe.', filename_asc);
        continue;
      end

    else
      fprintf('ASC-Datei %s existiert bereits, überspringe Konvertierung.\n', filename_asc);
    end


    disp('read asc data into FieldTrip ''raw'' data')
    cfg = [];
    cfg.dataset          = filename_asc;
    cfg.montage.tra      = eye(7);
    cfg.montage.labelorg = {'1', '2', '3', '4', '5', '6', '7'}; % 7 chans
    cfg.montage.labelnew = {'EYE_TIMESTAMP', 'EYE_L_HORIZONTAL', 'EYE_L_VERTICAL', 'EYE_L_DIAMETER', 'EYE_R_HORIZONTAL', 'EYE_R_VERTICAL', 'EYE_R_DIAMETER'};
    data = ft_preprocessing(cfg);

    filename_mat = [eyename '.mat'];
    disp(filename_mat)
    save(filename_mat, 'data')

    plotit = 0;
    if ispc && plotit
      if isfield(data, 'trial') && ~isempty(data.trial)
        cfg = [];
        cfg.channel = [4 7];
        ft_databrowser(cfg, data);
      else
        warning('data leer oder ungültig → kein ft_databrowser');
      end

      figure; plot(data.trial{1}(2,:), data.trial{1}(3,:)); hold on
      plot(data.trial{1}(5,:), data.trial{1}(6,:))
    end

    % Events einlesen
    disp 'read events from asc'
    event = ft_read_event(filename_asc);
    event_ori = event;
    event = struct2table(event); % convert to table
    event.type = string(event.type);
    event.value = string(event.value);

    if ispc && plotit
      if isfield(data_L, 'trial') && ~isempty(data_L.trial)
        cfg = [];
        cfg.viewmode       = 'vertical';
        cfg.channel        = { 'EYE_L_DIAMETER' 'EYE_R_DIAMETER' };
        cfg.preproc.demean = 'yes';
        cfg.event          = event_ori(strcmp({event_ori.type}, 'BLINK'));
        ft_databrowser(cfg, data_L);
      else
        warning('data_L leer oder ungültig → ft_databrowser wird übersprungen.');
      end
    end

    % BlinkArtefakte interpolieren
    disp('interpolate blinks separately for L and R eye')
    cfg = [];
    cfg.channel = { 'EYE_L_DIAMETER'};
    data_L = ft_selectdata(cfg, data);
    cfg.channel = {'EYE_R_DIAMETER'};
    data_R = ft_selectdata(cfg, data);
    cfg.channel = {'EYE_TIMESTAMP' 'EYE_L_HORIZONTAL', 'EYE_L_VERTICAL', 'EYE_R_VERTICAL', 'EYE_R_HORIZONTAL' };
    data_timestamp = ft_selectdata(cfg, data);

    % put EyeLink blink messages into FieldTrip compatible format
    inds = event.type == 'BLINK' & event.value == 'L';
    padpre = 150; % 150 samples = 0.15 s with 1000 Hz sampling rate
    padpost = 200; % 150 samples = 0.15 s with 1000 Hz sampling rate
    movement_L = [event.sample(inds)-padpre event.sample(inds)+event.duration(inds)+padpost];
    movement_L = movement_L( movement_L(:,1) > 1,:); % remove blink if already blinking at run start
    movement_L = movement_L(movement_L(:,2) <= event.sample(end),:);
    inds = event.type == 'BLINK' & event.value == 'R';
    movement_R = [event.sample(inds)-padpre event.sample(inds)+event.duration(inds)+padpost];
    movement_R = movement_R( movement_R(:,1) > 1,:);
    movement_R = movement_R(movement_R(:,2) <= event.sample(end),:);

    % plot the artifacts in databrowser
    if ispc && plotit
      cfg = [];
      cfg.artfctdef.blinks.artifact = movement_L;
      cfg.channel =  'EYE_L_DIAMETER';
      cfg.preproc.demean = 'yes'; % this makes the data zero-centered
      cfg = ft_databrowser(cfg, data_L);
    end


    % replace the artifacts with nans
    cfg=[];
    cfg.artfctdef.reject = 'nan';
    cfg.artfctdef.jump.artifact = movement_L;
    data_L = ft_rejectartifact(cfg, data_L);
    cfg.artfctdef.jump.artifact = movement_R;
    data_R = ft_rejectartifact(cfg, data_R);

    % append channels again
    data = ft_appenddata([], data_timestamp, data_L, data_R);

    event2 = event;
    event2(event2.value == "L",:) = [];
    event2(event2.value == "R",:) = [];

    if ispc && plotit
      cfg = [];
      cfg.event = table2struct(event2);
      cfg.channel =  'EYE_L_DIAMETER';
      cfg.preproc.demean = 'yes'; % this makes the data zero-centered
      ft_databrowser(cfg, data_L);
    end
    
    % append channels again
    data = ft_appenddata([], data_timestamp, data_L, data_R);
    
    event2 = event;
    event2(event2.value == "L",:) = [];
    event2(event2.value == "R",:) = [];
    
      if ispc && plotit
        cfg = [];
        cfg.event = table2struct(event2);
        cfg.channel =  'EYE_L_DIAMETER';
        cfg.preproc.demean = 'yes'; % this makes the data zero-centered
        ft_databrowser(cfg, data_L);
      end
    
    % plot the artifacts in databrowser
    if ispc && plotit
      if isfield(data, 'trial') && ~isempty(data.trial)
        cfg = [];
        cfg.channel        = {'EYE_L_DIAMETER' 'EYE_R_DIAMETER'} ;
        cfg.preproc.demean = 'yes';
        cfg.viewmode       = 'butterfly';
        cfg = ft_databrowser(cfg, data);
      else
        warning('data leer oder ungültig → kein butterfly-Plot.');
      end
    end

    % Assume your data is in FieldTrip format, e.g., data_pupil
    cfg = [];
    cfg.artfctdef.jump.channel = {'EYE_L_DIAMETER' 'EYE_R_DIAMETER'} ;  % Adjust if needed
    cfg.artfctdef.jump.medianfilter = 'yes';    % Smoother estimation (optional)
    cfg.artfctdef.jump.medianfiltord = 9;       % Window size for median filter
    cfg.artfctdef.jump.absdiff = 'yes';         % Use absolute difference (default)
    cfg.artfctdef.jump.cutoff = 10;             % Threshold for jump size (adjust!)
    cfg.artfctdef.jump.trlpadding = 0;
    cfg.artfctdef.jump.artpadding = 0.1;
    cfg.artfctdef.jump.fltpadding = 0;
    if plotit
      cfg.artfctdef.jump.interactive = 'yes';     % Optional GUI
    end
    % Run artifact detection
    [~, artifact_jump] = ft_artifact_jump(cfg, data);


    % replace the artifacts with nans
    cfg=[];
    cfg.artfctdef.reject = 'nan';
    cfg.artfctdef.jump.artifact = artifact_jump;
    data = ft_rejectartifact(cfg, data);

    cfg=[];
    cfg.channel = {'EYE_L_DIAMETER' 'EYE_R_DIAMETER'} ;
    data = ft_selectdata(cfg, data);

    % Make sure that data does not start or end with nan, before
    % interpolation
    data = ft_fill_leading_trailing_nans(data);

    % linear interpolation of nan timepoints
    cfg=[];
    cfg.prewindow = 1./data.fsample;
    cfg.postwindow = 1./data.fsample;
    % cfg.method      ='pchip';
    datakeep = data;
    data = ft_interpolatenan(cfg, data);
    % plot the artifacts in databrowser
    if ispc && plotit
      if isfield(data, 'trial') && ~isempty(data.trial)
        cfg = [];
        cfg.channel        = {'EYE_L_DIAMETER' 'EYE_R_DIAMETER'} ;
        cfg.preproc.demean = 'yes';
        cfg.viewmode       = 'butterfly';
        cfg = ft_databrowser(cfg, data);
      else
        warning('data leer oder ungültig → kein butterfly-Plot.');
      end
    end

    % plot the artifacts in databrowser
    if ispc && plotit
      cfg = [];
      % cfg.artfctdef.blinks.artifact = artifact_threshold;
      cfg.channel =  'EYE_L_DIAMETER';
      cfg.preproc.demean = 'yes'; % this makes the data zero-centered
      cfg = ft_databrowser(cfg, data);
    end

    % define trials
    disp('define trials')
    cfg=[];
    cfg.trialdef.trg = 'stim';
    cfg.trialdef.begtim = -5;
    cfg.trialdef.endtim = 10;
    cfg.event = event;
    cfg.fsample = data.fsample;
    cfg.trialfun = 'sortTrials';
    %cfg.trialdata = Data.Task(3).Results{1};
    cfg = ft_definetrial(cfg); % make trl matrix
    data = ft_redefinetrial(cfg, data); %make trials

    % Bedingung zuweisen anhand Dateiname
    if contains(lower(filename_asc), 'konservativ')
      condcode = 1;
      cond = "conservative";
    elseif contains(lower(filename_asc), 'liberal')
      condcode = 2;
      cond = "liberal";
    elseif contains(lower(filename_asc), 'baseline')
      condcode = 0;
      cond = "baseline";
    elseif contains(lower(filename_asc), 'training')
      condcode = 3;
      cond = "training";
    else
      condcode = "unknown";
    end


    nTrials = numel(data.trial);
    % data.trialinfo = table(Data.Task(3).Results{1}, 'Correct', 'Stimulus'); % add trial information
    % data.trialinfo.condition = 'Liberal'; % figure out which condition, e.g. based on filename, contains('conservative')
    % data.trialinfo.condition = 'Conservative';
    data.trialinfo = table();
    data.trialinfo.condition = repmat(cond, nTrials,1);
    data.trialinfo.runno = repmat(irun, nTrials,1);
    % data.trialinfo = [data.trialinfo, repmat(condcode, nTrials, 1)];
    % 
    % % ggf. Zeitverlauf berechnen
    % timelock = ft_timelockanalysis([], data);
    % % cfg = [];
    % % cfg.keeptrials = 'yes';
    % % timelock = ft_timelockanalysis(cfg, data); % make timelock: timelock.trial matrix with dims rpt_chan matrix, easier to handle
    % all_timelock{end+1} = timelock;


    if ispc && plotit
      cfg=[]; cfg.channel = {'EYE_L_DIAMETER' 'EYE_R_DIAMETER'};
      ft_singleplotER(cfg, timelock)
    end

    alldata{end+1} = data;
    
  end % irun
  
  data_avg{icond} = ft_appenddata([], alldata{:});
  timelock_avg{icond} = ft_timelockanalysis([], data);

  cfg = [];
  cfg.resample = 'yes';
  cfg.resamplefs = 100;
  cfg.detrend = 'no';
  data_avg{icond} = ft_resampledata(cfg,  data_avg{icond});
  timelock_avg{icond} = ft_resampledata(cfg,  timelock_avg{icond});

  outputfile = string(SubjectCode) + "_" + conds(icond) + ".mat";
  data = data_avg{icond};
  timelock = timelock_avg{icond};
  save(fullfile(outputFolder, outputfile), 'data', 'timelock');

end
if plotit
  ft_singleplotER([], timelock_avg{:})
  legend(conds)
end


% %% Nur konservative + liberale Runs zusammenfassen
% data_cl = {};
% for i = 1:numel(alldata)
%   thisData = alldata{i};
%   if isfield(thisData, 'trialinfo') && any(thisData.trialinfo(:,end) == 1 | thisData.trialinfo(:,end) == 2)
%     data_cl{end+1} = thisData;
%   end
% end
% 
% if ~isempty(data_cl)
%   cfg = [];
%   data_cl_merged = ft_appenddata(cfg, data_cl{:});
%   outfile_cons_lib = [SubjectCode '_allruns.mat'];
%   fprintf('Speichere konservativ+liberal zusammengefasste Datei: %s\n', outfile_cons_lib);
%   save(outfile_cons_lib, 'data_cl_merged');
% else
%   warning('Keine konservativen/liberalen Daten zum Zusammenfassen gefunden.');
% end
% 
% % plotting
% cons = all_timelock{3};
% lib = all_timelock{7};
% 
% % if ispc && plotit
% figure
% cfg=[];
% cfg.channel = {'EYE_L_DIAMETER' 'EYE_R_DIAMETER'};
% % cfg.baseline = [-1 0];
% % cfg.baselinetype = 'relchange';
% ft_singleplotER(cfg, cons, lib) %
% xline(0)
% legend({'cons' 'lib'})
% xlabel('Time from stim onset')
% ylabel('Pupil raw')
% if size(cons.avg, 1) >= 4
%   figure; plot(cons.time, cons.avg(4,:)); xline(0); xticks(-4:0.25:10); xline(1.01)
% else
%   warning('Es gibt nur %d Kanal(e); Kanal 4 ist nicht vorhanden.', size(cons.avg, 1));
% end
% 
% % end
% 
% cfg=[];
% cfg.keepsampleinfo = 'no';
% data = ft_appenddata(cfg, alldata{:});
% clear alldata
% 
% disp(outfile)
% save(outfile, '-struct', 'data')

