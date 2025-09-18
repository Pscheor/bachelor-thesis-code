function [trl, event] = sortTrials(cfg)
% make trl, columns:  
%  1) start 
%  2) end
%  3) offset 
%  4) study phase (1) / test phase (2)
%  5) stimulus category: fractals (1) / landscapes (2) / naturals1 (3) / streets1 (4) /  streets2 (5)
%  6) Pic number of image presented
%  7) quizpic number (nan for test)
%  8) quizpic: old (2) /new (1)
%  9) response: yes old (2) / yes new (1)
% 10) response button: left (1) / right (2)
% 11) RT for quizpic response
% 12) trigger wrt scanning onset ('fMRI scanning Start time')
% 13) C1 scale 1-8 mean % HMAX values (load in vals from file)
% 14) C2 patch 1-8 mean % HMAX values (load in vals from file)
% 15) runno in chronological order
% 16) remembered at test yes (1, Hit) or no (0, Miss)
% 17) RT at test
% 18) trial counter per run (1:30) 

% hdr    = cfg.headerfile;
fsr    = cfg.fsample;         % in Hz
% trg    = cfg.trialdef.trg;    % 'stim' or 'resp' or 'baseline'
begtrl = cfg.trialdef.begtim; % in seconds
endtrl = cfg.trialdef.endtim; % in seconds
event = cfg.event;
% runinfo = cfg.runinfo;
% trialdata = cfg.trialdata;

trialind = event.value == "Sound: no " | event.value == "Sound: on ";
% trialind = event.value == "Sound: on ";
% trialind = event.value == "Sound: no ";
% trialind = contains(event.value, "Response") & not(contains(event.value, "Start"));

% trialind = contains(event.value, "Confidence");
% trialind = contains(event.value, "BaselineBeforeRun");
trgval = event.value(trialind);
trgsmp = event.sample(trialind);
trgtimestamp = event.timestamp(trialind);

% valsplit = split(trgval);
% RT = str2double(valsplit(:,3));
% RT = RT(RT < 2.5);
% mean(RT) % 1.0192

% 1st timestamp: 2646778


ntrials = sum(trialind);
trl = zeros(ntrials, 3);
trl(:,1) = round(trgsmp + begtrl*fsr);  
trl(:,2) = round(trgsmp + endtrl*fsr); 
trl(:,3) = round(begtrl*fsr);

% trl = trl(2:end,:);

% trl = [trl trialdata];
% 
% trgval = {event( find(strcmp('Trial',{event.type})) ).value};
% trgsmp = [event( find(strcmp('Trial',{event.type})) ).sample];
% 
% abortedrun = 0;
% if runinfo.subjno == 10 && runinfo.runno == 5
%   trgval = trgval(1:144); % end of file exception, recovered with -failsafe
%   trgsmp = trgsmp(1:144);
%   abortedrun = 1;
% end
% 
% disp(runinfo.hand)
% 
% % find stim onset triggers
% stimstartend = strfind(trgval, ['Stim ' runinfo.category]); 
% stimstartend_index = find(not(cellfun('isempty', stimstartend)));
% stimstart = strfind(trgval, 'Start');
% stimstart_index = find(not(cellfun('isempty', stimstart)));
% trig_ind_stim = intersect(stimstartend_index,stimstart_index);
% 
% % create trl matrix
% trl = zeros(length(trig_ind_stim), 11);
% trl(:,1) = round(trgsmp(trig_ind_stim)' + begtrl*fsr); % find stim onset samples
% trl(:,2) = round(trgsmp(trig_ind_stim)' + endtrl*fsr); % trial end 
% trl(:,3) = round(begtrl*fsr);
% 
% exp_phases = {'study' 'test'};
% trl(:,4) = find( strcmp(exp_phases, runinfo.phase ));
% 
% % TODO fix for naturals2 subjects: only take first naturals
% stimcategories = {'fractals' 'landscapes' 'naturals1' 'streets1' 'streets2'};
% trl(:,5) = find(strcmp(stimcategories, runinfo.category ));
% 
% % get stim pic nr from trgval strings
% if trl(1,5) == 2 || trl(1,5) == 5 % for landscapes and streets2 format is imageXXX 
%     pat = [exp_phases{trl(1,4)} ' image(\w*) Start'];
% else
%     pat = [exp_phases{trl(1,4)} ' (\w*) Start']; 
% end
% temp = regexp(trgval(trig_ind_stim)', pat, 'tokens'); 
% temp = [temp{:}];   
% temp = [temp{:}];
% trl(:,6) = cellfun(@str2num, temp)';
% 
% % get resp triggers
% trig_ind_resp = find(cellfun(@isempty, strfind(trgval, 'Response')) == 0);
% % keep only first resp index if subject responded twice in one trial
% for i = 2:length(trig_ind_resp)
%     if trig_ind_resp(i) == trig_ind_resp(i-1)+1
%         trig_ind_resp(i) = NaN;
%     end
% end
% for i = 3:length(trig_ind_resp)
%     if trig_ind_resp(i) == trig_ind_resp(i-2)+2
%         trig_ind_resp(i) = NaN;
%     end
% end
% for i = 4:length(trig_ind_resp)
%     if trig_ind_resp(i) == trig_ind_resp(i-3)+3
%         trig_ind_resp(i) = NaN;
%     end
% end
% for i = 5:length(trig_ind_resp)
%     if trig_ind_resp(i) == trig_ind_resp(i-4)+4
%         trig_ind_resp(i) = NaN;
%     end
% end
% for i = 6:length(trig_ind_resp)
%     if trig_ind_resp(i) == trig_ind_resp(i-5)+5
%         trig_ind_resp(i) = NaN;
%     end
% end
% 
% nan_ind = find(isnan(trig_ind_resp));
% for i = 1:length(nan_ind)
%     trig_ind_resp = [trig_ind_resp(1:nan_ind(i)-1), trig_ind_resp(nan_ind(i)+1:end)];
%     nan_ind = nan_ind-1;
% end
% if isnan(trig_ind_resp(end))
%     trig_ind_resp = trig_ind_resp(1:end-1);
% end
% 
% % get quiz indices
% trig_ind_quiz = find(cellfun(@isempty,strfind(trgval, 'Stim quiz')) == 0);
% % if numel(trig_ind_quiz) ~= 30 % in case trigger missing
% %   trig_ind_quiz = find(cellfun(@isempty,strfind(trgval, 'Stim Quiz')) == 0);
% % end
% 
% % find missing trial indices
% trig_ind_miss = trig_ind_quiz + 1;
% trig_ind_miss = find(ismember(trig_ind_miss, trig_ind_resp) == 0);
% 
% trig_ind_resp_complete = trig_ind_resp;
% for i = 1:length(trig_ind_miss)
%     trig_ind_resp_complete = [trig_ind_resp_complete(1:trig_ind_miss(i)-1), NaN, trig_ind_resp_complete(trig_ind_miss(i):end)];
% end
% 
% 
% % % get quiz pic nr from trgval strings
% % if trl(1,5) == 2 || trl(1,5) == 5 % for landscapes and streets2 there is imageXXX
% %     pat = [stimcategories{trl(1,5)} ' image(\w*) RT'];
% % else
% %     pat = [stimcategories{trl(1,5)} ' (\w*) RT'];
% % end
% % temp = regexp(trgval(trig_ind_resp)', pat, 'tokens'); 
% % temp = [temp{:}];
% % temp = [temp{:}];
% 
% 
% % alternative way of finding quiz pic from trgval strings
% % above does not work if responses are missing, 
% % because trials without response are not in trgval
% 
% 
% if trl(1,5) == 2 || trl(1,5) == 5 % for landscapes and streets2 there is imageXXX
%     pat = [' image(\w*) Start']; % filter changed to 'imagexxx Start'
% else
%     pat = [' (\w*) Start'];
% end
% 
% temp = regexp(trgval(trig_ind_quiz)', pat, 'tokens');
% % if isempty([temp{:}])
% %   if trl(1,5) == 2 || trl(1,5) == 5 % for landscapes and streets2 there is imageXXX
% %     pat = [' image(\w*) End']; % filter changed to 'imagexxx Start'
% %   else
% %     pat = [' (\w*) End'];
% %   end
% %   temp = regexp(trgval(trig_ind_quiz)', pat, 'tokens');
% % end
% temp = [temp{:}];
% temp = [temp{:}];
% 
% 
% if trl(1,4) == 1 % study phase
%     trl(:,7) = cellfun(@str2num, temp)'; 
% else
%     trl(:,7) = NaN; % no quiz pic in test phase
%     pic_resps = cellfun(@str2num, temp)'; % but check for omissions (max 10 s RT for test)    
%     omissions = ismember(trl(:,6), pic_resps); % 0 for omissions; is the value in trl column 6 found somewhere in pics_resp
% end
% 
% %Get correctness
% if trl(1,4) == 2 % test phase
%     studytrialinfo = cfg.studytrialinfo(cfg.studytrialinfo(:,2) == trl(1,5),:); % select category from study trialinfo
%     correct_ans = ismember(trl(:,6), studytrialinfo(:,3)) + 1; % NEW=1, OLD=2
% else
%     correct_ans = (trl(:,6) == trl(:,7)) + 1; 
% end
% % trl(:,9) = trl(:,8) == correct_ans; % 0 is incorrect, 1 = correct TODO put old/new in col 8
% trl(:,8) = correct_ans; % put actual old/new in col 8
% 
% 
% % get response and response button
% if trl(1,4) == 1 % study phase
%     pat = 'Response (\w*) Quiz'; 
% else
%     pat = 'Response (\w*) Test';    % for getting resp button out 
% end
% temp = regexp(trgval(trig_ind_resp)', pat, 'tokens'); 
% temp = [temp{:}];   
% temp = [temp{:}];
% resps = zeros(length(trig_ind_resp),1); % for response
% button = zeros(length(trig_ind_resp),1); % for response button
% if runinfo.hand == 1
%     resps(find(strcmp('g', temp ))) = 1;
%     resps(find(strcmp('z', temp ))) = 2;
% else
%     resps(find(strcmp('z', temp ))) = 1;
%     resps(find(strcmp('g', temp ))) = 2;
% end
% button (find(strcmp('z', temp ))) = 1;
% button (find(strcmp('g', temp ))) = 2;
% 
% % align resps and button with trials where there was a response
% resps_complete = resps(1:length(trig_ind_resp));
% for i = 1:length(trig_ind_miss)
%     resps_complete = [resps_complete(1:trig_ind_miss(i)-1); NaN; resps_complete(trig_ind_miss(i):end)];
% end
% 
% button_complete = button(1:length(trig_ind_resp));
% for i = 1:length(trig_ind_miss)
%     button_complete = [button_complete(1:trig_ind_miss(i)-1); NaN; button_complete(trig_ind_miss(i):end)];
% end
% 
% 
% %line up responses with corresponding stims
% if trl(1,4) == 2 % test phase
%     resps_linedup = zeros(length(trl),1);
%     resps_linedup(omissions,:) = resps(resps>0);
%     trl(:,9) = resps_linedup; 
% else
%     trl(:,9) = resps_complete; % no omissions possible for study
% end
% 
% if trl(1,4) == 2 % test phase
%     resps_linedup = zeros(length(trl),1);
%     resps_linedup(omissions,:) = resps(resps>0);
%     trl(:,10) = resps_linedup; 
% else
%     trl(:,10) = button_complete; % no omissions possible for study
% end
% 
% 
% % Get RT
% pat = 'RT (..\d*)'; 
% temp = regexp(trgval(trig_ind_resp)', pat, 'tokens'); % get RT from trgval strings
% temp = [temp{:}];   
% temp = [temp{:}];
% temp = round(cellfun(@str2num, temp)' * 1000);
% 
% % align RTs with trials where there was a response
% for i = 1:length(trig_ind_miss)
%     temp = [temp(1:trig_ind_miss(i)-1); NaN; temp(trig_ind_miss(i):end)];
% end
% 
% if trl(1,4) == 2 % test phase
%     trl(omissions,11) = temp;
% else
%     trl(:,11) = temp;
% 
% end
% 
% % stim onset triggers, same as trialonsets.txt
% % get fMRI scanning Start time
% scanstartind = cellfun(@(x) ~isempty(strfind(x, 'fMRI scanning Start time')), {event.value});
% scanstartsmp = event(scanstartind).sample;
% 
% trl(:,12) = trl(:,1)/1000 - scanstartsmp/1000 + 3 - 12 + 5; 
% % -3 for begtrl (see above), 12 s (TRs) removed from fMRI data BUT 'fMRI scanning Start time' only sent after 5 s
% % see EMcheckfMRItimings()
% 
% % load Hmax
% % also in: /Users/kloosterman/gridmaster2012/projectdata/eyemem/D_paradigm/stimuli_640x480/hmax
% if ismac
% %   load(fullfile('/Users/kloosterman/gridmaster2012/LNDG/EyeMem/study_information/D_paradigm/stimuli_640x480/hmax', ['hmax_' runinfo.category]))
%   load(fullfile('/Users/kloosterman/gridmaster2012/projectdata/eyemem/D_paradigm/stimuli_640x480/hmax', ['hmax_' runinfo.category]))
% else
%   load(fullfile('/home/mpib/LNDG/EyeMem/study_information/D_paradigm/stimuli_640x480/hmax', ['hmax_' runinfo.category]))
% end
% hmaxout.c2median = cell2mat(hmaxout.c2median);
% 
% [test,ind] = sort(trl(:,6)); % image1 was shown in trial 19
% % [trl(ind,6) hmaxout.picno]
% if abortedrun == 0
%   if corr(trl(ind,6), hmaxout.picno) ~= 1; error('Hmax pics do not match!'); end
%   % trl(ind,13) = hmaxout.c1median(:,1); % trl(ind,6) gives ascending images
%   trl(ind,13) = mean(hmaxout.c1median,2); % trl(ind,6) gives ascending images
%   trl(ind,14) = mean(hmaxout.c2median,2); % trl(ind,6) gives ascending images
% else
%   trl(ind,13) = mean(hmaxout.c1median(1:29),2); % trl(ind,6) gives ascending images
%   trl(ind,14) = mean(hmaxout.c2median(1:29),2); % trl(ind,6) gives ascending images
% end  
% trl(:,15) = runinfo.runno;
% % trl(ind,15) = hmaxout.c2median(:,1); 
% % trl(ind,16) = hmaxout.c2median(:,8); 
% % % [trl(ind,13)  hmaxout.c1median(:,1)]
% % % [trl(:,13)  hmaxout.c1median(:,1)]
% 
% % load behavior
% if ismac
%   load('/Users/kloosterman/gridmaster2012/projectdata/eyemem/preproc/behavior/Eyemem_behavior.mat');
% else
%   load('/mnt/beegfs/home/kloosterman/projectdata/eyemem/preproc/behavior/Eyemem_behavior.mat')
% end
% % behav_test = behav.test(runinfo.subjno).singletrial; %  condition target_present response accuracy RT picno
% behav_test = behavior.singletrial{runinfo.subjno}{2}; %  condition target_present response accuracy RT picno
% behav_test = behav_test(behav_test(:,1) == trl(1,5),:);
% 
% trlind = ismember(behav_test(:,6) , trl(:,6)); %
% behav_test = behav_test(trlind,:); % select trials from test shown at study
% 
% % sort behav_test trials to match order of study
% for itrial = 1:size(trl,1)
%   testdat = behav_test(behav_test(:,6) == trl(itrial,6),:);  
%   if isempty(testdat) 
%     warning('pic not presented during test?? TODO find out')
%     continue;
%   end
%   trl(itrial,16) = testdat(4); % add col remembered (1) or not (0)
%   trl(itrial,17) = testdat(5);% add col rt at test
%   % TODO add col running nr
% %   trl(itrial,18) = find(behav_test(:,6) == trl(itrial,6));
% end
% 
% trl(:,18) = 1:size(trl,1); % to check matching eye and fMRI
% 
